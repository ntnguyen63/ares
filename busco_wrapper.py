from typing import List, Optional, Union, Dict

import subprocess

from pathlib import Path


class BuscoResult:
    def __init__(
        self,
        *,
        assembly: str,
        busco_path: Optional[str] = None,
        busco_score: Optional[float] = None,
        lineage: Optional[str] = None,
    ):
        self.assembly = assembly
        self.busco_score = busco_score
        self.busco_path = busco_path
        self.lineage = lineage


def initialize_busco_result(assembly: Union[str, BuscoResult]) -> BuscoResult:
    # start as a busco result from previous polishing
    if isinstance(assembly, BuscoResult):
        return assembly

    # or as a unpolished assembly
    return BuscoResult(
        assembly=assembly,
    )


def is_improved(*, new_busco: BuscoResult, old_busco: BuscoResult):
    if new_busco.busco_score > old_busco.busco_score:
        return True
    else:
        return False


def is_first(busco_result: BuscoResult):
    if None in (busco_result.busco_score, busco_result.busco_path):
        return True
    else:
        return False


def get_busco_score(short_summary):
    """get busco Complete score from short_summary.txt"""

    with open(short_summary) as f:
        for line in f:
            line = line.strip()
            if line.startswith(("#", "*")) or line == "":
                continue
            elif line.startswith("C:"):
                line = line.replace("%", "").replace("[", ",").replace("]", "")
                return float(line.split(",")[0].split(":")[1])


def run_busco(assembly: str, outdir: str, lineage: str):
    outdir = Path(outdir)

    # no slashes allowed in -o parameter so put stem as output
    stem = outdir.stem
    subprocess.run(
        f"busco -m genome -i {assembly} -o {stem} -l {lineage} --cpu 30",
        shell=True,
    )

    if "/" in outdir.as_posix():
        subprocess.run(f"mv {stem} {outdir}", shell=True)

    path = list(outdir.glob("*/short_summary.txt"))
    if not path:
        raise Exception("cannot find short_summary.txt")

    short_summary = path[0]
    busco_score = get_busco_score(short_summary)
    return BuscoResult(
        assembly=assembly,
        busco_score=busco_score,
        busco_path=outdir,
        lineage=lineage,
    )
