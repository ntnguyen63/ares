import os
import shutil
import subprocess

from pathlib import Path


def test_polishing():
    p = subprocess.Popen(
        "python ares.py test/drafs test/long.fasta test/r1.fastq test/r2.fastq".split()
    )
    stdout, stderr = p.communicate()
    assert p.returncode == 0
    assert os.path.exists("best/polish.fasta")


import polisher


def test_pilon_polish():
    pipeline = polisher.PolishPipeline(
        root_dir="polishing",
        drafts=[Path("test/drafts/assembly.fasta")],
        long_reads="test/long.fasta",
        short_reads=["test/r1.fasta", "test/r2.fasta"],
        threads=25,
        lineage="fungi_odb10",
    )
    Path("polishing/assembly").mkdir(exist_ok=True, parents=True)
    busco_result = pipeline.pilon_polish(
        "polishing/assembly", "test/drafts/assembly.fasta"
    )
    assert busco_result.assembly
    assert busco_result.busco_score > 30
    assert busco_result.busco_path


def test_is_improved():
    better = polisher.BuscoResult(
        assembly="test/drafts/assembly.fasta", busco_score=89, busco_path=None
    )
    worst = polisher.BuscoResult(
        assembly="test/drafts/assembly.fasta", busco_score=70, busco_path=None
    )
    assert polisher.is_improved(new_busco=better, old_busco=worst)
    assert not polisher.is_improved(new_busco=worst, old_busco=better)


def test_is_first():
    busco_result = polisher.BuscoResult(
        assembly="test/drafts/assembly.fasta",
        busco_score=None,
        busco_path=None,
    )
    assert polisher.is_first(busco_result)


def test_run_busco():
    root_dir = "busco_test"
    Path(root_dir).mkdir(exist_ok=True)

    out = f"{root_dir}/busco_out"
    assembly = "test/drafts/assembly.fasta"
    lineage = "basidiomycota_odb10"
    if os.path.exists(out):
        shutil.rmtree(out, ignore_errors=True)

    result = polisher.run_busco(assembly, out, lineage)
    assert os.path.exists(out)
    assert result.assembly
    assert result.busco_score
    assert result.busco_path
