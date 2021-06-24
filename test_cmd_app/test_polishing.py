import os
import shutil
import subprocess

from pathlib import Path

from ..lib.polisher import (
    PolishPipeline,
    create_sorted_aln,
    pilon,
    ShortReadPolishRunner,
)

draft = Path("test_data/drafts/assembly.fasta")
r1, r2 = ["test_data/r1.fastq", "test_data/r2.fastq"]
long_reads = "test_data/long.fastq"
busco_lineage = "fungi_odb10"
threads = 30

pipeline = PolishPipeline(
    root_dir="polishing",
    drafts=[Path(draft)],
    long_reads=long_reads,
    short_reads=[r1, r2],
    threads=threads,
    lineage=busco_lineage,
)


def test_polish():
    Path("polishing").mkdir(exist_ok=True)
    busco_result = pipeline.polish(draft)
    assert os.path.exists(busco_result.assembly)
    assert busco_result.busco_score
    assert busco_result.busco_path


def test_medaka_polish():
    Path("polishing/assembly").mkdir(exist_ok=True, parents=True)
    busco_result = pipeline.medaka_polish(
        "polishing/assembly", draft.as_posix()
    )
    assert os.path.exists(busco_result.assembly)
    assert busco_result.busco_path
    assert busco_result.busco_score


def test_racon_polish():
    Path("polishing/assembly").mkdir(exist_ok=True, parents=True)
    busco_result = pipeline.racon_polish(
        "polishing/assembly", draft.as_posix()
    )
    assert os.path.exists(busco_result.assembly)
    assert busco_result.busco_path
    assert busco_result.busco_score


def test_create_sorted_aln():
    sorted_aln = create_sorted_aln(
        assembly=draft,
        r1=r1,
        r2=r2,
        threads=threads,
        out="sorted_aln",
    )
    assert sorted_aln
    assert sorted_aln.is_file()


def test_pilon_polish():
    Path("polishing/assembly").mkdir(exist_ok=True, parents=True)
    busco_result = pipeline.pilon_polish(
        "polishing/assembly", draft.as_posix()
    )
    assert busco_result.assembly
    assert busco_result.busco_score > 30
    assert busco_result.busco_path


def test_pilon():
    Path("polishing/assembly").mkdir(exist_ok=True, parents=True)
    busco_result = pilon(
        outdir="polishing/assembly",
        draft=draft.as_posix(),
        r1=r1,
        r2=r2,
        busco_lineage=busco_lineage,
        threads=threads,
    )
    assert busco_result.assembly
    assert busco_result.busco_score > 30
    assert busco_result.busco_path


def test_short_read_polish_runner():
    Path("polishing/assembly").mkdir(exist_ok=True, parents=True)
    runner = ShortReadPolishRunner(
        root_dir="polishing/assembly",
        draft=draft.as_posix(),
        r1=r1,
        r2=r2,
        lineage=busco_lineage,
        threads=threads,
        rounds=4,
    )
    best = runner.run(pilon)
    assert best.assembly
    assert best.busco_score
    assert best.busco_path
    assert os.path.exists("polishing/assembly/pilon")
