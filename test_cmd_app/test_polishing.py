import os
import shutil
import subprocess

from pathlib import Path

from ..polisher import PolishPipeline, create_sorted_aln

pipeline = PolishPipeline(
    root_dir="polishing",
    drafts=[Path("test_data/drafts/assembly.fasta")],
    long_reads="test_data/long.fastq",
    short_reads=["test_data/r1.fastq", "test_data/r2.fastq"],
    threads=25,
    lineage="fungi_odb10",
)


def test_polish():
    Path("polishing").mkdir(exist_ok=True)
    draft = Path("test_data/drafts/assembly.fasta")
    busco_result = pipeline.polish(draft)
    assert os.path.exists(busco_result.assembly)
    assert busco_result.busco_score
    assert busco_result.busco_path


def test_racon_polish():
    Path("polishing/assembly").mkdir(exist_ok=True, parents=True)
    busco_result = pipeline.racon_polish(
        "polishing/assembly", "test_data/drafts/assembly.fasta"
    )
    assert os.path.exists(busco_result.assembly)
    assert busco_result.busco_path
    assert busco_result.busco_score


def test_create_sorted_aln():
    sorted_aln = create_sorted_aln(
        assembly=Path("test_data/drafts/assembly.fasta"),
        r1="test_data/r1.fastq",
        r2="test_data/r2.fastq",
        threads=25,
        out="sorted_aln",
    )
    assert sorted_aln
    assert sorted_aln.is_file()


def test_pilon_polish():
    Path("polishing/assembly").mkdir(exist_ok=True, parents=True)
    busco_result = pipeline.pilon_polish(
        "polishing/assembly", "test_data/drafts/assembly.fasta"
    )
    assert busco_result.assembly
    assert busco_result.busco_score > 30
    assert busco_result.busco_path
