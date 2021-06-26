from typing import List, Optional, Union, Dict, Callable, Protocol

import subprocess

from pathlib import Path

from . import busco_wrapper

import time


def get_best_busco_result(results: List[busco_wrapper.BuscoResult]):
    return max(
        results,
        key=lambda result: result.busco_score,
    )


def create_dir_for_best(best: busco_wrapper.BuscoResult, outdir: str):
    subprocess.run(f"mkdir -p {outdir}", shell=True)
    subprocess.run(f"cp {best.assembly} {outdir}/polish.fasta", shell=True)
    subprocess.run(f"cp -r {best.busco_path} {outdir}/", shell=True)


def create_sorted_aln(
    *, assembly: Path, r1: str, r2: str, threads: int, out: str
) -> Path:
    subprocess.run(
        f"minimap2 -ax sr {assembly} {r1} {r2} | samtools view -u | samtools sort -@ {threads} > {out}",
        shell=True,
    )
    subprocess.run(f"samtools index {out}", shell=True)
    return Path(out)


class ShortReadPolisher(Protocol):
    def __call__(
        self,
        *,
        outdir: str,
        draft: str,
        r1: str,
        r2: str,
        lineage: str,
        threads: int,
    ) -> busco_wrapper.BuscoResult:
        ...


def pilon(
    *,
    outdir: str,
    draft: str,
    r1: str,
    r2: str,
    lineage: str,
    threads: int,
    rounds: int,
) -> busco_wrapper.BuscoResult:
    sorted_out = "sorted.bam"
    alignment = create_sorted_aln(
        assembly=Path(draft),
        r1=r1,
        r2=r2,
        threads=threads,
        out=sorted_out,
    )
    subprocess.run(
        [
            "pilon",
            "--fix",
            "all,amb",
            "--genome",
            draft,
            "--frags",
            alignment,
            "--threads",
            str(threads),
            "--outdir",
            outdir,
        ]
    )
    # some clean up
    subprocess.run(f"rm {sorted_out.split('.')[0]}.*", shell=True)
    polish = f"{outdir}/pilon.fasta"

    if not Path(polish).is_file():
        raise Exception(f"pilon was unable to polish {draft}")

    return busco_wrapper.run_busco(polish, f"{outdir}/busco_out", lineage)
    
def create_bam_bwa(*, genome: str, threads: int, read1: str, read2: str):
    subprocess.run(["bwa-mem2", "index", genome])
    subprocess.run( f"bwa-mem2 mem -t {threads} {genome} {read1} {read2}| \
                    samtools view --threads {threads} -F 0x4 -b -|samtools fixmate -m --threads {threads}  - -| \
                    samtools sort -m 8g --threads {threads} -|samtools markdup --threads {threads} -r - sorted.bam", shell=True
                   )
    subprocess.run(["samtools", "index", "sorted.bam"])
    subprocess.run(["samtools", "faidx", genome])
    return
    
def nextpolish(
    *,
    outdir: str,
    draft: str,
    r1: str,
    r2: str,
    lineage: str,
    threads: int,
    rounds: int,
) -> busco_wrapper.BuscoResult:
    sorted_out = "sorted.bam"
    create_bam_bwa(genome=draft, read1=r1, read2=r2, threads=threads)
    command = [
            "nextpolish1.py",
            "-g",
            draft,
            "-t",
            str(rounds+1), #does not accept 0
            "-p",
            str(threads),
            "-s",
            "sorted.bam", #it shows  PosixPath('sorted.bam') if using alignment arg?
            ">",
            outdir+'/nextpolish.fasta',
        ]
    print(command)
    Path(outdir.split('/')[0]+"/pilon/round_"+str(rounds)).mkdir(parents=True, exist_ok=True)
    #nextpolish has to be run in a shell to cat output to a file
    subprocess.run(f"nextpolish1.py -g {draft} -t 1 -p {str(threads)} -s sorted.bam > {outdir}/temp1.nextpolish.fasta", shell=True)
    # some clean up
    subprocess.run(f"rm {sorted_out.split('.')[0]}.*", shell=True)
    subprocess.run(f"rm {draft.split('/')[-1]}.*",shell=True)
    #redo polishing second mode
    create_bam_bwa(
        genome=outdir+'/temp1.nextpolish.fasta',
        read1=r1,
        read2=r2,
        threads=threads,
    )
    subprocess.run(f"nextpolish1.py -g {outdir+'/temp1.nextpolish.fasta'} -t 2 -p {str(threads)} -s sorted.bam > {outdir}/nextpolish.fasta", shell=True)
    # some clean up
    subprocess.run(f"rm {sorted_out.split('.')[0]}.*", shell=True)
    subprocess.run(f"rm {outdir}/temp1.nextpolish.fasta*", shell=True)
    polish = f"{outdir}/nextpolish.fasta"

    if not Path(polish).is_file():
        raise Exception(f"nextpolish was unable to polish {draft}")

    return busco_wrapper.run_busco(polish, f"{outdir}/busco_out", lineage)


class ShortReadPolishRunner:
    def __init__(
        self,
        root_dir: str,
        draft: str,
        r1: str,
        r2: str,
        threads: int,
        lineage: str,
        rounds: int,
    ):
        self.root_dir = root_dir
        self.draft = draft
        self.r1 = r1
        self.r2 = r2
        self.threads = threads
        self.lineage = lineage
        self.rounds = rounds
        self.all_busco_runs: List[busco_wrapper.BuscoResult] = []

    def run(self, polisher: ShortReadPolisher):
        busco_result = busco_wrapper.initialize_busco_result(self.draft)
        for i in range(self.rounds):
            pilon_out = f"{self.root_dir}/pilon/round_{i}"

            new_busco = polisher(
                outdir=pilon_out,
                draft=busco_result.assembly,
                r1=self.r1,
                r2=self.r2,
                threads=self.threads,
                lineage=self.lineage,
                rounds=i
            )

            self.all_busco_runs.append(new_busco)

            if busco_wrapper.is_first(
                busco_result
            ) or busco_wrapper.is_improved(
                new_busco=new_busco, old_busco=busco_result
            ):
                busco_result = new_busco
            else:
                break

        return busco_result


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
        results: List[busco_wrapper.BuscoResult] = [
            self.polish(draft) for draft in self.drafts if draft.is_file()
        ]
        best = max(
            results,
            key=lambda result: result.busco_score,
        )
        best_dir = f"{self.root_dir}/best"
        subprocess.run(f"mkdir -p {best_dir}", shell=True)
        subprocess.run(
            f"cp {best.assembly} {best_dir}/polish.fasta", shell=True
        )
        subprocess.run(f"cp -r {best.busco_path} {best_dir}/", shell=True)

        # create tsv file from all busco runs
        busco_wrapper.summarize_busco_runs(
            outdir=self.root_dir, best=best, busco_results=self.all_busco_runs
        )

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
            command = f"medaka_consensus -i {self.long_reads} -d {busco_result.assembly} -o {out_dir} -t {self.threads} -m r941_prom_fast_g303"
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
                f"pilon --fix all,amb --genome {busco_result.assembly} --frags {sorted_aln} --threads {self.threads} --outdir {pilon_out}",
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
