import os
import shutil
import subprocess

from pathlib import Path

from ..lib.busco_wrapper import (
    run_busco,
    is_improved,
    is_first,
    BuscoResult,
    summarize_busco_runs,
)

better = BuscoResult(
    assembly="test_data/drafts/assembly.fasta",
    busco_score=89,
    busco_path="test/busco/path/a1.fasta",
    lineage="fungi_odb10",
)
worst = BuscoResult(
    assembly="test_data/drafts/assembly.fasta",
    busco_score=70,
    busco_path="test/busco/path/a2.fasta",
    lineage="fungi_odb10",
)


def test_summarize_busco_runs():
    outdir = "summarize_busco_runs_test"
    Path(outdir).mkdir(exist_ok=True)
    summarize_busco_runs(
        outdir=outdir,
        best=better,
        busco_results=[better, worst],
    )
    assert os.path.exists(outdir)
    assert os.path.exists(f"{outdir}/results.tsv")
    assert os.path.exists(f"{outdir}/best/polish.fasta")


def test_is_improved():
    assert is_improved(new_busco=better, old_busco=worst)
    assert not is_improved(new_busco=worst, old_busco=better)


def test_is_first():
    busco_result = BuscoResult(
        assembly="test_data/drafts/assembly.fasta",
        busco_score=None,
        busco_path=None,
    )
    assert is_first(busco_result)


def test_run_busco():
    root_dir = "busco_test"
    Path(root_dir).mkdir(exist_ok=True)
    out = f"{root_dir}/busco_out"
    assembly = "test_data/drafts/assembly.fasta"
    lineage = "fungi_odb10"

    if os.path.exists(out):
        shutil.rmtree(out, ignore_errors=True)

    result = run_busco(assembly, out, lineage)
    assert os.path.exists(out)
    assert result.assembly
    assert result.busco_score
    assert result.busco_path
