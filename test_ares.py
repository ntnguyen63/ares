import os
import shutil
import subprocess


def test_polishing():
    p = subprocess.Popen(
        "python ares.py test/drafs test/long.fasta test/r1.fastq test/r2.fastq".split()
    )
    stdout, stderr = p.communicate()
    assert p.returncode == 0
    assert os.path.exists("best/polish.fasta")


import polisher


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
    assembly = "test/drafts/assembly.fasta"
    lineage = "basidiomycota_odb10"
    out = "busco_test"
    if os.path.exists(out):
        shutil.rmtree(out, ignore_errors=True)

    result = polisher.run_busco(assembly, out, lineage)
    assert os.path.exists(out)
    assert result.contigs
    assert result.busco_score
    assert result.busco_path
