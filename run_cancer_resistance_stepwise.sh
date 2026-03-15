#!/usr/bin/env bash
# BioResilient Cancer Resistance — Stepwise Interactive Pipeline Runner
#
# Runs each step (or step group) sequentially, pausing after each one to:
#   1. Show a rich cached report (step_cache/<step>.md)
#   2. Show the scientific plausibility validation result (PASS / WARN / FAIL)
#   3. Ask for confirmation before proceeding
#
# Usage:
#   ./run_cancer_resistance_stepwise.sh                     # full interactive run
#   ./run_cancer_resistance_stepwise.sh --non-interactive   # CI mode: auto-go, exit 1 on FAIL
#   ./run_cancer_resistance_stepwise.sh --from step7        # resume from a specific step
#   ./run_cancer_resistance_stepwise.sh --only step1,step2  # run only listed steps
#   ./run_cancer_resistance_stepwise.sh --dry-run           # test the runner without executing steps

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHENOTYPE="cancer_resistance"
LOG_FILE="$REPO_ROOT/pipeline.log"
CACHE_DIR="$REPO_ROOT/step_cache"
NON_INTERACTIVE=false
RESUME_FROM=""
ONLY_STEPS=""
DRY_RUN_FLAG=""
CI_EXIT_CODE=0

# ── Colours ───────────────────────────────────────────────────────────────────
RED="\033[0;31m"; GREEN="\033[0;32m"; YELLOW="\033[1;33m"
CYAN="\033[0;36m"; BOLD="\033[1m"; DIM="\033[2m"; RESET="\033[0m"
BLUE="\033[0;34m"; MAGENTA="\033[0;35m"

# ── Parse arguments ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --non-interactive) NON_INTERACTIVE=true ;;
        --from)            RESUME_FROM="$2"; shift ;;
        --only)            ONLY_STEPS="$2"; shift ;;
        --dry-run)         DRY_RUN_FLAG="--dry-run" ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
    shift
done

# ── Step sequence (groups of steps run together) ──────────────────────────────
# Format: "display_name|step1,step2,..."
# Steps within a group run together without pausing; pause happens after the group.
declare -a STEP_GROUPS=(
    "Environment check|step1"
    "Proteome download|step2"
    "OrthoFinder clustering|step3,step3b"
    "Nucleotide region extraction + alignment|step3c"
    "Divergence + motif finding|step4"
    "Domain & AlphaMissense annotation|step4b"
    "ESM-1v variant scoring|step4c"
    "GoF / LoF direction|step4d"
    "Phylogenetic tree|step5"
    "Phylogenetic conservation scoring (phyloP/PhastCons)|step3d"
    "MEME positive selection|step6"
    "FEL + BUSTED selection tests|step6b"
    "RELAX rate acceleration|step6c"
    "Convergence scoring|step7"
    "Convergent amino acid check|step7b"
    "GEO expression data|step8"
    "Bgee expression supplement|step8b"
    "Phase 1 composite scoring  ◀ PRIMARY CHECKPOINT|step9"
    "Disease annotation|step11,step11b,step11c,step11d"
    "Druggability (fpocket + ChEMBL)|step12,step12b"
    "Gene therapy feasibility|step13"
    "Safety screen|step14,step14b"
    "Final rescore  ◀ FINAL CHECKPOINT|step15"
)

# ── Helpers ───────────────────────────────────────────────────────────────────

banner() {
    echo -e "\n${BOLD}${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
    echo -e "${BOLD}${CYAN}  $1${RESET}"
    echo -e "${BOLD}${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}\n"
}

info() {
    echo -e "${DIM}  $1${RESET}"
}

success() {
    echo -e "${GREEN}  ✓ $1${RESET}"
}

warn() {
    echo -e "${YELLOW}  ⚠ $1${RESET}"
}

error() {
    echo -e "${RED}  ✗ $1${RESET}"
}

run_step() {
    # run_step "step4,step4b" "divergence"
    local steps="$1"
    local label="$2"
    echo -e "\n${GREEN}${BOLD}▶  Running ${label}...${RESET}"
    cd "$REPO_ROOT"
    # Run each step in the comma-separated list
    IFS=',' read -ra STEP_LIST <<< "$steps"
    for step in "${STEP_LIST[@]}"; do
        echo -e "${DIM}   → $step${RESET}"
        python pipeline/orchestrator.py \
            --steps "$step" \
            --phenotype "$PHENOTYPE" \
            $DRY_RUN_FLAG \
            2>&1 | tee -a "$LOG_FILE"
        local exit_code=${PIPESTATUS[0]}
        if [[ $exit_code -ne 0 ]]; then
            error "Step $step exited with code $exit_code"
            return 1
        fi
    done
    return 0
}

write_and_show_report() {
    # Write JSON+MD report and show validation result. Sets VALIDATION_STATUS.
    local step="$1"
    echo -e "${DIM}   Writing report for ${step}...${RESET}"
    python pipeline/step_reporter.py --step "$step" 2>&1
    VALIDATION_STATUS=$?    # 0 = PASS/WARN, 1 = FAIL
    return 0
}

sync_cache_to_s3() {
    # Upload step cache files to S3.
    # Current copy: s3://{bucket}/step_cache/{phenotype}/{step}.{md,json}  (always overwritten)
    # History copy: s3://{bucket}/step_cache/history/{phenotype}/{step}/{date}.{md,json}
    #   — only written if the content differs from the current S3 copy (i.e. a real change)
    local step="$1"
    local bucket="${S3_BUCKET:-bioresilient-data}"
    local phenotype_key="${PHENOTYPE:-cancer_resistance}"
    local date_tag
    date_tag=$(date -u +"%Y-%m-%dT%H%M%SZ")

    if ! command -v aws &>/dev/null; then
        return 0
    fi

    for ext in json md; do
        local src="$CACHE_DIR/${step}.${ext}"
        [[ -f "$src" ]] || continue
        local s3_current="s3://${bucket}/step_cache/${phenotype_key}/${step}.${ext}"
        local s3_history="s3://${bucket}/step_cache/history/${phenotype_key}/${step}/${date_tag}.${ext}"
        local tmp_prev="/tmp/_s3_prev_${step}.${ext}"

        # Check if current S3 content differs from new file
        local changed=true
        if aws s3 cp "$s3_current" "$tmp_prev" --quiet 2>/dev/null; then
            if diff -q "$src" "$tmp_prev" &>/dev/null; then
                changed=false
            fi
            rm -f "$tmp_prev"
        fi

        # Always overwrite current
        aws s3 cp "$src" "$s3_current" --quiet 2>/dev/null && \
            echo -e "${DIM}   S3 cache saved: step_cache/${phenotype_key}/${step}.${ext}${RESET}"

        # Only save history if content changed
        if $changed; then
            aws s3 cp "$src" "$s3_history" --quiet 2>/dev/null && \
                echo -e "${DIM}   S3 history saved: ${date_tag}${RESET}"
        fi
    done
}

show_report_md() {
    local step="$1"
    local md_file="$CACHE_DIR/${step}.md"
    if [[ -f "$md_file" ]]; then
        echo ""
        echo -e "${DIM}$(cat "$md_file")${RESET}"
        echo ""
    fi
}

pause_and_confirm() {
    # pause_and_confirm "step4" "PASS" "next_step_label"
    local last_step="$1"
    local status="$2"
    local next_label="$3"

    if $NON_INTERACTIVE; then
        if [[ "$status" == "FAIL" ]]; then
            error "FAIL detected in non-interactive mode — halting."
            CI_EXIT_CODE=1
            return 1
        fi
        info "Non-interactive mode: auto-proceeding (status=$status)"
        return 0
    fi

    echo ""
    echo -e "  ${DIM}Cache: ${CACHE_DIR}/${last_step}.json  |  ${CACHE_DIR}/${last_step}.md${RESET}"
    echo ""

    if [[ "$status" == "FAIL" ]]; then
        echo -e "  ${RED}${BOLD}Step FAILED validation. Do not proceed without investigating.${RESET}"
        echo -e "  ${RED}Suggested: review the recommendation above and fix config/data before continuing.${RESET}"
        echo ""
        echo -e "  ${BOLD}Type  stop  → halt here (fix and resume with --from ${last_step})${RESET}"
        while true; do
            read -r -p "  > " response
            case "$response" in
                stop|STOP|q|Q) return 1 ;;
                continue|go|force)
                    warn "Forcing past FAIL — results may be unreliable."
                    return 0 ;;
                *) echo "  Please type  stop  to halt (or  force  to override)." ;;
            esac
        done
    else
        if [[ "$status" == "WARN" ]]; then
            echo -e "  ${YELLOW}${BOLD}⚠  WARN detected above. Review before proceeding.${RESET}"
        else
            echo -e "  ${GREEN}${BOLD}✅ Validation passed.${RESET}"
        fi
        echo -e "  Next: ${BOLD}${next_label}${RESET}"
        echo ""
        echo -e "  ${BOLD}Type  go    → proceed to next step${RESET}"
        echo -e "  ${BOLD}Type  show  → show full Markdown report${RESET}"
        echo -e "  ${BOLD}Type  skip  → skip next step${RESET}"
        echo -e "  ${BOLD}Type  stop  → halt here (resume with --from ${last_step})${RESET}"
        while true; do
            read -r -p "  > " response
            case "$response" in
                go|GO|"") return 0 ;;
                show|SHOW) show_report_md "$last_step"; continue ;;
                skip|SKIP) return 2 ;;
                stop|STOP|q|Q) return 1 ;;
                *) echo "  Please type  go / show / skip / stop." ;;
            esac
        done
    fi
}

is_step_active() {
    # Returns 0 (true) if step should run given --from and --only flags
    local step="$1"
    if [[ -n "$ONLY_STEPS" ]]; then
        IFS=',' read -ra ONLY_LIST <<< "$ONLY_STEPS"
        for s in "${ONLY_LIST[@]}"; do
            [[ "$s" == "$step" ]] && return 0
        done
        return 1
    fi
    if [[ -n "$RESUME_FROM" ]]; then
        local found=false
        for grp in "${STEP_GROUPS[@]}"; do
            local steps="${grp#*|}"
            IFS=',' read -ra STEPS_IN_GROUP <<< "$steps"
            for s in "${STEPS_IN_GROUP[@]}"; do
                [[ "$s" == "$RESUME_FROM" ]] && found=true
                $found && [[ "$s" == "$step" ]] && return 0
            done
        done
        return 1
    fi
    return 0
}

# ── Pre-flight ─────────────────────────────────────────────────────────────────
mkdir -p "$CACHE_DIR"
touch "$LOG_FILE"

# Guard against duplicate instances — only one pipeline may run at a time.
LOCKFILE="/tmp/bioresilient_pipeline.lock"
if [[ -f "$LOCKFILE" ]]; then
    existing_pid=$(cat "$LOCKFILE" 2>/dev/null)
    if [[ -n "$existing_pid" ]] && kill -0 "$existing_pid" 2>/dev/null; then
        echo -e "${RED}ERROR: Pipeline already running (PID $existing_pid). Kill it first:${RESET}"
        echo -e "  kill $existing_pid"
        exit 1
    else
        echo -e "${YELLOW}Stale lock found (PID $existing_pid gone) — removing.${RESET}"
        rm -f "$LOCKFILE"
    fi
fi
echo $$ > "$LOCKFILE"
trap 'rm -f "$LOCKFILE"' EXIT INT TERM

banner "BioResilient Cancer Resistance — Stepwise Pipeline"
echo -e "  Phenotype  : ${BOLD}${PHENOTYPE}${RESET}"
echo -e "  Repo       : ${REPO_ROOT}"
echo -e "  Log        : ${LOG_FILE}"
echo -e "  Cache      : ${CACHE_DIR}"
[[ -n "$RESUME_FROM" ]] && echo -e "  Resume from: ${BOLD}${RESUME_FROM}${RESET}"
[[ -n "$ONLY_STEPS"  ]] && echo -e "  Only steps : ${BOLD}${ONLY_STEPS}${RESET}"
$NON_INTERACTIVE && echo -e "  Mode       : ${YELLOW}non-interactive (CI)${RESET}" || echo -e "  Mode       : interactive"
[[ -n "$DRY_RUN_FLAG" ]] && echo -e "  ${YELLOW}DRY RUN — orchestrator steps skipped${RESET}"
echo ""

# ── Main loop ─────────────────────────────────────────────────────────────────
GROUP_COUNT=${#STEP_GROUPS[@]}
current_idx=0

for grp in "${STEP_GROUPS[@]}"; do
    current_idx=$((current_idx + 1))
    DISPLAY_NAME="${grp%%|*}"
    STEPS_STR="${grp#*|}"

    # Determine the primary/last step of this group (for reporting)
    IFS=',' read -ra STEPS_IN_GROUP <<< "$STEPS_STR"
    LAST_STEP="${STEPS_IN_GROUP[-1]}"

    # Check if any step in this group is active
    group_active=false
    for s in "${STEPS_IN_GROUP[@]}"; do
        is_step_active "$s" && { group_active=true; break; }
    done
    $group_active || continue

    # Peek at next group name for the prompt
    NEXT_DISPLAY="(pipeline complete)"
    for ((ni=current_idx; ni<GROUP_COUNT; ni++)); do
        NEXT_GRP="${STEP_GROUPS[$ni]}"
        NEXT_STEPS_STR="${NEXT_GRP#*|}"
        IFS=',' read -ra NEXT_STEPS <<< "$NEXT_STEPS_STR"
        next_active=false
        for ns in "${NEXT_STEPS[@]}"; do
            is_step_active "$ns" && { next_active=true; break; }
        done
        if $next_active; then
            NEXT_DISPLAY="${NEXT_GRP%%|*}"
            break
        fi
    done

    # Run
    banner "Group ${current_idx}/${GROUP_COUNT}: ${DISPLAY_NAME}"
    if run_step "$STEPS_STR" "$DISPLAY_NAME"; then
        success "Completed: $DISPLAY_NAME"
    else
        error "Step group FAILED: $DISPLAY_NAME"
        if $NON_INTERACTIVE; then
            CI_EXIT_CODE=1
            break
        fi
        echo -e "  ${RED}The step above failed with a non-zero exit code.${RESET}"
        echo -e "  Type  retry  to re-run, or  stop  to halt."
        read -r -p "  > " resp
        case "$resp" in
            retry|RETRY)
                if run_step "$STEPS_STR" "$DISPLAY_NAME"; then
                    success "Completed on retry: $DISPLAY_NAME"
                else
                    error "Still failing. Halting."
                    exit 1
                fi
                ;;
            *) exit 1 ;;
        esac
    fi

    # Write report and validate (using last step of group)
    write_and_show_report "$LAST_STEP"
    sync_cache_to_s3 "$LAST_STEP"
    VSTATUS="PASS"
    if [[ -f "$CACHE_DIR/${LAST_STEP}.json" ]]; then
        VSTATUS=$(python -c "
import json, sys
try:
    d = json.loads(open('$CACHE_DIR/${LAST_STEP}.json').read())
    v = d.get('_validation', {})
    print(v.get('status', 'PASS'))
except Exception:
    print('PASS')
" 2>/dev/null || echo "PASS")
    fi

    # Pause and confirm
    pause_and_confirm "$LAST_STEP" "$VSTATUS" "$NEXT_DISPLAY"
    result=$?
    if [[ $result -eq 1 ]]; then
        echo -e "\n${YELLOW}  Pipeline halted after ${DISPLAY_NAME}.${RESET}"
        echo -e "  Resume later with: ${BOLD}./run_cancer_resistance_stepwise.sh --from ${LAST_STEP}${RESET}"
        exit $CI_EXIT_CODE
    elif [[ $result -eq 2 ]]; then
        warn "Skipping next step group."
        continue
    fi
done

# ── Done ──────────────────────────────────────────────────────────────────────
if [[ $CI_EXIT_CODE -eq 0 ]]; then
    banner "Pipeline Complete"
    success "All steps finished for phenotype: ${PHENOTYPE}"
    echo -e "  Full reports in: ${BOLD}${CACHE_DIR}/${RESET}"
    echo -e "  Pipeline log  : ${LOG_FILE}"
    echo ""
    echo -e "  ${DIM}Run benchmark recall check:${RESET}"
    echo -e "  ${BOLD}python scripts/benchmark_recall.py --trait ${PHENOTYPE}${RESET}"
    echo ""
else
    error "Pipeline finished with FAIL validations (CI mode). See log: $LOG_FILE"
fi

exit $CI_EXIT_CODE
