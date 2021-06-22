from typing import List, Optional, Union, Dict

import subprocess

from pathlib import Path

from . import busco_wrapper


class PolishPipeline:
    def __init__(
        self,
        *,
        root_dir: str,
        drafts: List[Path],
        long_reads: str,
        short_reads: List[str],
        threads: int,
        lineage: str,
    ):
        self.root_dir = root_dir
        self.drafts = drafts
        self.long_reads = long_reads
        self.r1, self.r2 = short_reads
        self.MEDAKA_ROUNDS = 4
        self.PILON_ROUNDS = 4
        self.threads = threads
        self.busco_lineage = lineage

        self.all_busco_runs: List[busco_wrapper.BuscoResult] = []

    def run(self):
        busco_results: List[busco_wrapper.BuscoResult] = [
            self.polish(draft) for draft in self.drafts if draft.is_file()
        ]
        best_polish = max(
            busco_results, key=lambda busco_result: busco_result.busco_score
        )
        best_dir = f"{self.root_dir}/best"
        subprocess.run(f"mkdir -p {best_dir}", shell=True)
        subprocess.run(
            f"cp {best_polish.assembly} {best_dir}/polish.fasta", shell=True
        )

        # create tsv file with results
        with open(f"{self.root_dir}/results.tsv", "w") as results:
            headers: List[str] = [
                "assembly",
                "lineage",
                "busco_score",
                "busco_path",
                "is_best",
            ]
            headers_as_str: str = "\t".join(headers)
            results.write(f"{headers_as_str}\n")
            for busco_result in self.all_busco_runs:
                row = "\t".join(
                    [
                        str(busco_result.__dict__[header])
                        for header in headers
                        if header != "is_best"
                    ]
                )
                row.append(
                    True
                    if busco_result.assembly == best_polish.assembly
                    else False
                )
                results.write(f"{row}\n")

    def polish(self, draft: Path):
        polish_dir = f"{self.root_dir}/{draft.stem}"
        Path(polish_dir).mkdir(exist_ok=True)
        medaka_polish = self.medaka_polish(polish_dir, draft)
        return self.pilon_polish(polish_dir, draft)

    def medaka_polish(self, polish_dir, draft) -> busco_wrapper.BuscoResult:
        busco_result = self.racon_polish(polish_dir, draft)
        self.all_busco_runs.append(busco_result)

        for i in range(self.MEDAKA_ROUNDS):
            out_dir = f"{polish_dir}/medaka/round_{i}"
            command = f"medaka_consensus -i {self.long_reads} -d {busco_result.assembly} -o {out_dir} -t {self.threads} -m r941_min_fast_g303"
            subprocess.run(command, shell=True)

            polish = f"{out_dir}/consensus.fasta"
            if not Path(polish).is_file():
                raise Exception(
                    f"medaka was unable to polish {busco_result.assembly} on round {i}"
                )

            new_busco_result = busco_wrapper.run_busco(
                polish, f"{out_dir}/busco_out", self.busco_lineage
            )

            self.all_busco_runs.append(new_busco_result)

            # just ran first polish or the nth polish has improved assembly
            if busco_wrapper.is_first(
                busco_result
            ) or busco_wrapper.is_improved(
                new_busco=new_busco_result, old_busco=busco_result
            ):
                busco_result = new_busco_result
            else:
                break

        return busco_result

    def pilon_polish(self, polish_dir, draft) -> busco_wrapper.BuscoResult:
        busco_result = busco_wrapper.initialize_busco_result(draft)
        self.all_busco_runs.append(busco_result)

        for i in range(self.PILON_ROUNDS):
            sorted_out = "sorted.bam"
            pilon_out = f"{polish_dir}/pilon/round_{i}"
            sorted_aln = create_sorted_aln(
                assembly=Path(busco_result.assembly),
                r1=self.r1,
                r2=self.r2,
                threads=self.threads,
                out=sorted_out,
            )
            subprocess.run(
                f"pilon --genome {busco_result.assembly} --frags {sorted_aln} --threads {self.threads} --outdir {pilon_out}",
                shell=True,
            )
            # some clean up
            subprocess.run(f"rm {sorted_out.split('.')[0]}.*", shell=True)
            polish = f"{pilon_out}/pilon.fasta"

            if not Path(polish).is_file():
                raise Exception(
                    f"pilon was unable to polish {busco_result.assembly} on round {i}"
                )

            new_busco_result = busco_wrapper.run_busco(
                polish, f"{pilon_out}/busco_out", self.busco_lineage
            )
            self.all_busco_runs.append(new_busco_result)

            # just ran busco for the best time
            if busco_wrapper.is_first(
                busco_result
            ) or busco_wrapper.is_improved(
                new_busco=new_busco_result, old_busco=busco_result
            ):
                busco_result = new_busco_result
            else:
                break

        return busco_result

    def racon_polish(self, polish_dir, draft) -> busco_wrapper.BuscoResult:
        busco_result = busco_wrapper.initialize_busco_result(draft)

        racon_out = Path(f"{polish_dir}/racon")
        racon_out.mkdir(exist_ok=True)
        racon_polish = Path(f"{racon_out}/racon.fasta")
        paf = "minimap2.racon.paf"
        subprocess.run(
            f"minimap2 {busco_result.assembly} {self.long_reads} > {paf}",
            shell=True,
        )

        subprocess.run(
            f"racon -t {self.threads} -m 8 -x -6 -g -8 -w 500 {self.long_reads} {paf} {busco_result.assembly} > {racon_polish}",
            shell=True,
        )
        if not racon_polish.is_file():
            raise Exception(
                f"Unable to polish {busco_result.assembly} with racon"
            )
        return busco_wrapper.run_busco(
            racon_polish.as_posix(),
            f"{racon_out}/busco_out",
            self.busco_lineage,
        )


def create_sorted_aln(
    *, assembly: Path, r1: str, r2: str, threads: int, out: str
) -> Path:
    subprocess.run(
        f"minimap2 -ax sr {assembly} {r1} {r2} | samtools view -u | samtools sort -@ {threads} > {out}",
        shell=True,
    )
    subprocess.run(f"samtools index {out}", shell=True)
    return Path(out)
