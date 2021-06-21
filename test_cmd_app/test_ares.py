import os
import shutil
import subprocess

from pathlib import Path


def test_areas():
    p = subprocess.Popen(
        "python ares.py ../test_data/drafs ../test_data/long.fasta ../test_data/r1.fastq ../test_data/r2.fastq".split()
    )
    stdout, stderr = p.communicate()
    assert p.returncode == 0
    assert os.path.exists("best/polish.fasta")
