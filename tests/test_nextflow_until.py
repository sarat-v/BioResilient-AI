from pathlib import Path
import subprocess


REPO_ROOT = Path(__file__).resolve().parents[1]


def run_nextflow_preview(*args: str, dag_path: Path) -> str:
    command = [
        "nextflow",
        "run",
        "nextflow/main.nf",
        "-profile",
        "local",
        "-preview",
        "-with-dag",
        str(dag_path),
        *args,
    ]
    result = subprocess.run(
        command,
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout + result.stderr


def test_until_step4d_stops_after_phase1_sequence(tmp_path):
    dag_path = tmp_path / "until-step4d.dot"
    run_nextflow_preview("--until", "step4d", dag_path=dag_path)

    dag = dag_path.read_text()

    assert "PHASE1_SEQUENCE:variant_direction" in dag
    assert "PHASE1_EVOLUTION:build_species_tree" not in dag
    assert "PHASE1_EXPRESSION:expression_analysis" not in dag
    assert "PHASE2_CLINICAL:alphagenome_regulatory" not in dag
