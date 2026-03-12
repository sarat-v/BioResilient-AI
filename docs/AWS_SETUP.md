# BioResilient AWS Setup Guide

This guide walks you through setting up the cloud infrastructure for one complete
18-species cancer resistance run. Total estimated cost: **~$9–11** on Spot pricing.

---

## What you need before starting

### Already done (you confirmed these)
- ✅ AWS account
- ✅ AWS CLI installed (`aws --version` should print a version)

### Check these now — on your **local Mac terminal**

**1. AWS CLI is configured with credentials:**
```bash
aws sts get-caller-identity
```
Expected output: your Account ID, UserId, and ARN. If it says "Unable to locate
credentials", run `aws configure` and enter your Access Key ID, Secret Access Key,
region (e.g. `us-east-1`), and output format (`json`).

**2. You have an EC2 key pair for SSH:**
```bash
aws ec2 describe-key-pairs --query 'KeyPairs[*].KeyName'
```
If empty, create one in **AWS Console → EC2 → Key Pairs → Create key pair**.
Download the `.pem` file to `~/.ssh/` and run:
```bash
chmod 400 ~/.ssh/<your-key-name>.pem
```

**3. Your AWS CLI region is set:**
```bash
aws configure get region
```
Should print something like `us-east-1`. All resources (VPC, EC2, RDS, S3) should
be in the same region to avoid cross-region data transfer fees.

---

## Where commands run — quick reference

| Section | Run commands on |
|---------|----------------|
| §1 Networking | AWS Console (browser) — no CLI commands needed |
| §2 RDS | AWS Console (browser) |
| §3 S3 | Your **local Mac terminal** (uses `aws` CLI) |
| §4 IAM | Your **local Mac terminal** (uses `aws` CLI), or Console |
| §5 Launch EC2 | AWS Console (browser) |
| §6 SSH + clone | Your **local Mac terminal** first, then the **EC2 terminal** |
| §7 Env vars | **EC2 terminal** |
| §8 Bootstrap | **EC2 terminal** |
| §9–11 Run pipeline | **EC2 terminal** |

---

## Architecture overview

```
Your Mac
  │
  │  (SSH)
  ▼
┌──────────────────────────────────────────────────────────────┐
│  EC2 Spot Instance (c6i.4xlarge)                             │
│  16 vCPU / 32 GB RAM / Ubuntu 22.04                         │
│  ┌──────────────────────────────────────────────┐           │
│  │  bioresillient conda env                      │           │
│  │  run_cancer_resistance_stepwise.sh            │    ───►  │  S3 bucket
│  │  step_cache/  (JSON + MD reports)             │           │  bioresillient-data
│  └──────────────────────────────────────────────┘           │  ├── proteomes/
│                      │                                       │  ├── genomes/
│                      │                                       │  ├── nucleotide_regions/
│                      │                                       │  ├── cache/
│                      ▼                                       │  └── step_cache/
│  RDS PostgreSQL 15 (db.t3.medium)                           └──────────────────
│  bioresillient database (private — only EC2 can reach it)
└──────────────────────────────────────────────────────────────┘
```

RDS is **not publicly accessible** — your Mac cannot connect to it directly.
The EC2 instance is the only thing that talks to RDS.

---

## Instance recommendations

### Phase 1 — Steps 1–9 (OrthoFinder, MAFFT, IQ-TREE2, HyPhy)

| Property | Value |
|----------|-------|
| Instance | **c6i.4xlarge** |
| vCPU / RAM | 16 vCPU / 32 GB RAM |
| On-demand | ~$0.68/hr |
| Spot price | ~$0.20–0.27/hr |
| Why | All Phase 1 tools are CPU-bound and multi-threaded. 32 GB handles 18-species OrthoFinder (~22 GB peak RAM). c6i (Ice Lake) gives best price/CPU ratio in AWS. |
| Est. time | ~16–20 h on Spot |
| Est. cost | **~$4–6 on Spot** |

### Step 4c only — ESM-1v scoring (GPU-bound, optional)

| Property | Value |
|----------|-------|
| Instance | g4dn.xlarge |
| vCPU / RAM | 4 vCPU / 16 GB RAM / T4 GPU 16 GB VRAM |
| Spot price | ~$0.16/hr |
| Why | ESM-1v 650M model needs ~6 GB VRAM; T4 fits safely. CPU inference takes ~8 h; T4 takes ~25 min. |
| Est. cost | ~$0.04 on Spot (25 min) |
| Option | Skip the instance switch — run step4c on c6i CPU. Adds ~8 h but saves the switch complexity. |

### Practical simplification

Run everything on a single **c6i.4xlarge** Spot instance, step 4c on CPU.
Total: **~$7–9** for one complete 18-species run.

---

## Required AWS services

| Service | Spec | Cost |
|---------|------|------|
| EC2 Spot | c6i.4xlarge, Ubuntu 22.04 LTS, 100 GB gp3 root | ~$6.50/run |
| RDS PostgreSQL 15 | db.t3.medium, 20 GB gp2, no Multi-AZ | ~$1.63/run |
| S3 | Private bucket, same region | ~$1.15 (persists for reruns) |
| IAM | EC2 instance role | free |

---

## Step-by-step setup

---

### 1. Networking
**Where: AWS Console (browser) → search "VPC" in the top search bar**

You need a VPC with a public subnet so your EC2 instance gets a public IP
(so you can SSH in), and a private subnet for RDS (so only EC2 can reach it).

**In AWS Console → VPC → Your VPCs → Create VPC:**
```
Resources to create:  VPC and more
Name tag:             bioresillient
IPv4 CIDR:            10.0.0.0/16
Number of AZs:        1
Public subnets:       1
Private subnets:      0   (RDS will go in the public subnet but blocked by SG)
NAT gateways:         None
```
Click **Create VPC**. This automatically creates the subnet, internet gateway, and route table.

**Create security groups — in AWS Console → VPC → Security Groups → Create:**

Security group 1 (for your EC2 instance):
```
Name:         bioresillient-ec2-sg
VPC:          bioresillient-vpc
Inbound rule: SSH | TCP | Port 22 | My IP   ← click "My IP", it fills automatically
Outbound:     All traffic (default)
```

Security group 2 (for RDS — no inbound from internet, only from EC2):
```
Name:         bioresillient-rds-sg
VPC:          bioresillient-vpc
Inbound rule: PostgreSQL | TCP | Port 5432 | Custom → bioresillient-ec2-sg
Outbound:     All traffic (default)
```

> **Why two SGs?** The RDS group only allows connections that come from the EC2
> security group. Your Mac, or anyone on the internet, cannot reach the database at all.

---

### 2. RDS — PostgreSQL 15
**Where: AWS Console → search "RDS" → Create database**

```
Engine:              PostgreSQL
Engine version:      PostgreSQL 15.x (latest 15)
Template:            Dev/Test
DB instance ID:      bioresillient
Master username:     bioresillient
Master password:     <create a strong password — save this, you'll need it later>
Instance class:      db.t3.medium
Storage:             20 GB gp2
Storage autoscaling: Off
VPC:                 bioresillient-vpc
Subnet group:        default (auto-created)
Public access:       No
VPC security group:  bioresillient-rds-sg  (remove the default sg)
Multi-AZ:            No
DB name:             bioresillient
```

Click **Create database**. It takes ~5 minutes to become available.

Once available, go to the database → **Connectivity & security** tab and note the
**Endpoint** — it looks like:
```
bioresillient.abc123xyz.us-east-1.rds.amazonaws.com
```
Save this. You will set it as `RDS_HOST` later.

---

### 3. S3 bucket
**Where: your local Mac terminal**

```bash
# Create the bucket (replace us-east-1 with your region if different)
aws s3 mb s3://bioresillient-data --region us-east-1

# Block all public access
aws s3api put-public-access-block \
    --bucket bioresillient-data \
    --public-access-block-configuration \
    "BlockPublicAcls=true,IgnorePublicAcls=true,BlockPublicPolicy=true,RestrictPublicBuckets=true"
```

Verify it was created:
```bash
aws s3 ls | grep bioresillient
```

> S3 bucket names are globally unique. If `bioresillient-data` is taken, use
> `bioresillient-data-<your-initials>` and update `S3_BUCKET` accordingly throughout.

---

### 4. IAM role for EC2
**Where: your local Mac terminal**

This role lets the EC2 instance read/write S3 and connect to RDS without needing
hard-coded credentials on the instance.

```bash
# Create the role
aws iam create-role \
    --role-name bioresillient-ec2-role \
    --assume-role-policy-document \
    '{"Version":"2012-10-17","Statement":[{"Effect":"Allow","Principal":{"Service":"ec2.amazonaws.com"},"Action":"sts:AssumeRole"}]}'

# Attach S3 and RDS permissions
aws iam attach-role-policy \
    --role-name bioresillient-ec2-role \
    --policy-arn arn:aws:iam::aws:policy/AmazonS3FullAccess

aws iam attach-role-policy \
    --role-name bioresillient-ec2-role \
    --policy-arn arn:aws:iam::aws:policy/AmazonRDSFullAccess

# Create an instance profile and attach the role
aws iam create-instance-profile \
    --instance-profile-name bioresillient-ec2-profile

aws iam add-role-to-instance-profile \
    --instance-profile-name bioresillient-ec2-profile \
    --role-name bioresillient-ec2-role
```

Verify:
```bash
aws iam get-instance-profile --instance-profile-name bioresillient-ec2-profile \
    --query 'InstanceProfile.Roles[0].RoleName'
# Should print: "bioresillient-ec2-role"
```

---

### 5. Launch EC2 Spot instance
**Where: AWS Console → EC2 → Launch instances**

```
Name:             bioresillient-phase1
AMI:              Ubuntu Server 22.04 LTS (64-bit x86)
                  (search "ubuntu 22.04" in the AMI search — use the one from Canonical)
Instance type:    c6i.4xlarge
Key pair:         <select your existing key pair>
Network:          bioresillient-vpc
Subnet:           the public subnet created in §1
Auto-assign IP:   Enable
Security group:   bioresillient-ec2-sg
Storage:          100 GB gp3  (change from the default 8 GB)
IAM profile:      bioresillient-ec2-profile
```

**To use Spot pricing (saves ~60%):**
Expand **Advanced details** → scroll to **Purchasing option** → check
**Request Spot instances** → set max price to `0.68` (on-demand price as ceiling).

Click **Launch instance**.

Once running (takes ~1 min), go to the instance → copy the **Public IPv4 address**.

---

### 6. SSH in and clone the repo
**Where: your local Mac terminal first, then you switch to the EC2 terminal**

**On your Mac:**
```bash
ssh -i ~/.ssh/<your-key-name>.pem ubuntu@<ec2-public-ip>
```

You are now inside the EC2 instance. Everything from here runs **on EC2**.

**On EC2:**
```bash
# Update packages
sudo apt-get update && sudo apt-get install -y git

# Clone the repo
git clone https://github.com/<your-org>/bioresillient.git
cd bioresillient
```

> **Tip:** Keep this SSH session open. If you disconnect, use
> `ssh -i ~/.ssh/<key>.pem ubuntu@<ip>` again and `cd bioresillient` to resume.
> For long runs it's worth using `tmux` or `screen` so the pipeline keeps running
> even if your laptop closes:
> ```bash
> sudo apt-get install -y tmux
> tmux new -s pipeline
> # (to reattach later: tmux attach -t pipeline)
> ```

---

### 7. Set environment variables
**Where: EC2 terminal**

Open `~/.bashrc` in a text editor:
```bash
nano ~/.bashrc
```

Add these lines at the bottom (fill in your actual values):
```bash
export DEPLOYMENT=cloud
export RDS_HOST=bioresillient.xxxx.us-east-1.rds.amazonaws.com   # from §2
export RDS_PASSWORD=<the password you set in §2>
export S3_BUCKET=bioresillient-data
export NCBI_API_KEY=<your-ncbi-key>    # free at ncbi.nlm.nih.gov/account — takes 30 sec
export NCBI_EMAIL=<your-email>
# Optional — pipeline works without these:
export ALPHAGENOME_API_KEY=<your-key>
export ANTHROPIC_API_KEY=<your-key>
```

Save (`Ctrl+O`, `Enter`, `Ctrl+X`) then reload:
```bash
source ~/.bashrc
```

Verify the key ones are set:
```bash
echo $RDS_HOST
echo $S3_BUCKET
```

---

### 8. Bootstrap environment
**Where: EC2 terminal, inside the `bioresillient` repo directory**

```bash
bash scripts/setup_cloud.sh
```

This single script does everything:
1. Installs Miniconda
2. Creates the `bioresillient` conda environment from `environment.yml`
3. Installs all bioinformatics tools via conda/bioconda
4. Writes `config/environment.yml` from your env vars
5. Runs Alembic DB migrations against RDS
6. Seeds the species registry into the DB
7. Creates all S3 directory prefixes
8. Validates all tools are present

This takes **~15–20 minutes** on first run (conda solve + large package downloads).

**Tools installed:**

| Tool | Used in | Notes |
|------|---------|-------|
| OrthoFinder | step 3 | ortholog clustering |
| DIAMOND | step 3 | OrthoFinder backend |
| MAFFT | step 4 | protein MSA |
| IQ-TREE2 | step 5 | phylogenetic tree |
| HyPhy | steps 6–6c | selection tests |
| fpocket | step 12 | druggability |
| minimap2 | step 3c | nucleotide alignment |
| LASTZ | step 3c | fallback aligner |
| PHAST (phyloP, phastCons) | step 3d | phylogenetic conservation scoring |

`minimap2`, `LASTZ`, and `PHAST` are optional — the pipeline degrades gracefully
without them (steps 3c/3d produce NULL scores, validation emits `INFO` not `FAIL`).

---

### 9. Activate conda and verify
**Where: EC2 terminal**

The setup script activates conda automatically, but if you start a new SSH session:
```bash
conda activate bioresillient
```

Verify the DB connection works:
```bash
python -c "
from db.session import get_session
with get_session() as s:
    print('DB OK:', s.execute('SELECT 1').scalar())
"
```
Should print: `DB OK: 1`

Verify S3 access:
```bash
aws s3 ls s3://bioresillient-data/
```
Should list the prefixes created by the setup script (`proteomes/`, `genomes/`, etc.).

---

### 10. Verify the full environment (step1 gate)
**Where: EC2 terminal**

```bash
python pipeline/orchestrator.py --steps step1 --phenotype cancer_resistance
```

Expected output: all tools found ✓, DB ping < 100ms, no FAIL gates.

If any tool shows NOT FOUND, install it manually:
```bash
conda install -c bioconda <toolname>
```

---

### 11. Start the pipeline
**Where: EC2 terminal (run inside tmux if you haven't already)**

```bash
./run_cancer_resistance_stepwise.sh
```

The pipeline pauses after every step and shows a prompt:
```
  ✅ Validation passed.
  Next: OrthoFinder clustering

  Type  go    → proceed to next step
  Type  show  → show full Markdown report
  Type  stop  → halt here (resume with --from step2)
  >
```

Type `go` to proceed. Type `stop` to pause and review before continuing.

**To resume after stopping or a Spot interruption:**
```bash
./run_cancer_resistance_stepwise.sh --from <last-completed-step>
```
e.g. `--from step5` — the pipeline picks up from there, skipping already-done steps.

---

## Cost summary (single c6i.4xlarge Spot run)

| Resource | Duration | Rate | Cost |
|----------|----------|------|------|
| c6i.4xlarge Spot | ~24 h | $0.27/hr | **~$6.50** |
| RDS db.t3.medium | ~24 h | $0.068/hr | **~$1.63** |
| S3 storage | ~50 GB | $0.023/GB | **~$1.15** |
| Data transfer | minimal | | **~$0.10** |
| **Total** | | | **~$9.38** |

S3 storage persists between runs — reruns only pay for compute (~$8.13).

> **To avoid surprise charges:** stop the EC2 instance when not running the pipeline
> (AWS Console → EC2 → select instance → Instance state → Stop). RDS is billed by
> the hour even when idle — consider stopping it too if you won't run for several days
> (RDS → Databases → select → Actions → Stop temporarily).

---

## Resuming after a Spot interruption

Spot instances can be interrupted with a 2-minute warning. The pipeline is
designed to survive this cleanly:

1. S3 holds all proteomes, aligned orthogroups, and the treefile permanently
2. RDS holds all intermediate results (orthologs, motifs, scores)
3. Relaunch a new Spot instance with the **same IAM profile and env vars** and run:

```bash
conda activate bioresillient
cd bioresillient
./run_cancer_resistance_stepwise.sh --from <last-completed-step>
```

The runner skips already-completed steps and picks up where it left off.

---

## Instance-to-instance handoff for GPU step

If you want to run step 4c on a GPU instance to save ~8 h of CPU time:

```bash
# On c6i instance — run up to step 4b then stop
./run_cancer_resistance_stepwise.sh --only "step1,step2,step3,step3b,step3c,step3d,step4,step4b"

# Launch a g4dn.xlarge with the same IAM role and env vars
# On g4dn instance:
conda activate bioresillient && cd bioresillient
./run_cancer_resistance_stepwise.sh --only "step4c"

# Return to c6i (or continue on g4dn) for the rest
./run_cancer_resistance_stepwise.sh --from step4d
```

---

## Monitoring

All commands below run on **EC2 terminal**:

```bash
# Watch pipeline output in real time
tail -f pipeline.log

# See step reports as they are written (one .md file per step)
watch -n 30 ls -lt step_cache/

# Read the last completed step report
cat step_cache/<stepname>.md

# Check the database is growing as expected
psql -h $RDS_HOST -U bioresillient -d bioresillient \
    -c "SELECT schemaname, relname, pg_size_pretty(pg_total_relation_size(relid))
        FROM pg_catalog.pg_statio_user_tables
        ORDER BY pg_total_relation_size(relid) DESC LIMIT 10;"

# Check S3 data sizes
aws s3 ls s3://bioresillient-data/ --recursive --human-readable --summarize | tail -3
```
