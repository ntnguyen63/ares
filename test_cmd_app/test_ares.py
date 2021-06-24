import os
import shutil
import subprocess

from pathlib import Path


def test_ares():
    p = subprocess.Popen(
        "python ares.py test_data/drafts test_data/long.fastq test_data/r1.fastq test_data/r2.fastq fungi_odb10 30".split()
    )
    stdout, stderr = p.communicate()
    assert p.returncode == 0
    assert os.path.exists("ares_out")
    assert os.path.exists("ares_out/best/polish.fasta")


def test_pilon():
    p = subprocess.Popen(
        "python ares.py pilon pilon_command_test test_data/drafts/assembly.fasta test_data/r1.fastq test_data/r2.fastq fungi_odb10 30".split()
    )
    stdout, stderr = p.communicate()
    assert p.returncode == 0
    assert os.path.exists("pilon_command_test")
    assert os.path.exists("pilon_command_test/best/polish.fasta")
