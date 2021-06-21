from typing import List, Optional, Union

import subprocess

from pathlib import Path


class BuscoResult:
    def __init__(
        self,
        *,
        assembly: str,
        busco_path: Optional[str] = None,
        busco_score: Optional[float] = None,
    ):
        self.assembly = assembly
        self.busco_score = busco_score
        self.busco_path = busco_path


def initialize_busco_result(assembly: Union[str, BuscoResult]):
    # start polished draft from other polishing tool
    if isinstance(assembly, BuscoResult):
        return assembly

    # or new draft
    else:
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
    stdout, stderr = subprocess.run(
        f"busco -m genome -i {assembly} -o {stem} -l {lineage} --cpu 30",
        shell=True,
    )

    path = list(outdir.glob("*/short_summary.txt"))
    if not path:
        raise Exception("cannot find short_summary.txt")

    short_summary = path[0]
    busco_score = get_busco_score(short_summary)
    return BuscoResult(
        assembly=assembly, busco_score=busco_score, busco_path=outdir
    )


# class PolishPipeline:
#    def __init__(
#        self, *, root_dir, drafts, long_reads, short_reads, threads, lineage
#    ):
#        self.root_dir = root_dir
#        self.drafts = drafts
#        self.long_reads = long_reads
#        self.r1, self.r2 = short_reads
#        self.MEDAKA_ROUNDS = 4
#        self.PILON_ROUNDS = 4
#        self.threads = threads
#        self.busco_lineage = lineage
#
#    def run(self):
#        busco_results = [self.polish(draft) for draft in self.drafts]
#        best_polish = max(
#            busco_results, key=lambda busco_result: busco_result.busco_score
#        )
#        best_dir = f"{self.root_dir}/best"
#        shell("mkdir -p {best_dir}")
#        shell("cp {best_polish.contigs} {best_dir}/polish.fasta")
#
#    def polish(self, draft):
#        polish_dir = f"{self.root_dir}/{draft.stem}"
#        Path(polish_dir).mkdir(exist_ok=True)
#
#        # pilon first
#        pilon_polish = self.pilon_polish(polish_dir, draft)
#
#        # then medaka
#        medaka_polish = self.medaka_polish(polish_dir, pilon_polish.contigs)
#
#        return medaka_polish
#
#    def medaka_polish(self, polish_dir, draft):
#        busco_result = initialize_busco_result(draft)
#
#        for i in range(self.MEDAKA_ROUNDS):
#            out_dir = f"{polish_dir}/medaka/round_{i}"
#            command = f"medaka_consensus -i {self.long_reads} -d {busco_result.contigs} -o {out_dir} -t {self.threads} -m r941_min_fast_g303"
#            subprocess.run(command.split(" "))
#
#            polish = f"{out_dir}/consensus.fasta"
#            if not Path(polish).is_file():
#                raise Exception(
#                    f"medaka was unable to polish {busco_result.contigs} on round {i}"
#                )
#
#            new_busco_result = run_busco(
#                polish, f"{out_dir}/busco_out", self.busco_lineage
#            )
#
#            # just ran first polish or the nth polish has improved assembly
#            if is_first(busco_result) or is_improved(
#                new_busco_result, busco_result
#            ):
#                busco_result = new_busco_result
#            else:
#                break
#
#            return busco_result
#
#    def pilon_polish(self, polish_dir, draft):
#        busco_result = initialize_busco_result(draft)
#
#        for i in range(self.PILON_ROUNDS):
#            sorted_aln = "sorted.bam"
#            pilon_out = f"{polish_dir}/pilon/round_{i}"
#            shell(
#                "minimap2 -ax sr {busco_result.contigs} {self.r1} {self.r2} | samtools view -u | samtools sort -@ {self.threads} > {sorted_aln}"
#            )
#            shell("samtools index {sorted_aln}")
#            shell(
#                "pilon --genome {busco_result.contigs} --frags {sorted_aln} --threads {self.threads} --outdir {pilon_out}"
#            )
#            shell("rm {sorted_aln}")
#            polish = f"{pilon_out}/pilon.fasta"
#
#            if not Path(polish).is_file():
#                raise Exception(
#                    f"pilon was unable to polish {busco_result.contigs} on round {i}"
#                )
#
#            new_busco_result = run_busco(
#                polish, f"{pilon_out}/busco_out", self.busco_lineage
#            )
#
#            # just ran busco for the best time
#            if is_first(busco_result) or is_improved(
#                new_busco_result, busco_result
#            ):
#                busco_result = new_busco_result
#            else:
#                break
#
#        return busco_result
