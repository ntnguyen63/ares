import os
import shutil
import subprocess

from pathlib import Path

from ..busco_wrapper import run_busco, is_improved, is_first, BuscoResult


def test_is_improved():
    better = BuscoResult(
        assembly="test_data/drafts/assembly.fasta",
        busco_score=89,
        busco_path=None,
    )
    worst = BuscoResult(
        assembly="test_data/drafts/assembly.fasta",
        busco_score=70,
        busco_path=None,
    )
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
