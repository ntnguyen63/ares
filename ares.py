#!/usr/bin/env python

from pathlib import Path

import typer

from lib.polisher import (
    PolishPipeline,
    ShortReadPolishRunner,
    pilon,
)

from lib.busco_wrapper import summarize_busco_runs

app = typer.Typer()


def handle_drafts(dir: str):
    drafts_dir = Path(dir)
    if not drafts_dir.is_dir():
        raise Exception(f"{dir} is not a directory")

    drafts_filtered = [
        draft for draft in drafts_dir.glob("*.fasta") if draft.is_file()
    ]
    if not drafts_filtered:
        raise Exception(f"No fasta files in {dir}")
    return drafts_filtered


@app.command()
def ares(
    drafts: str, nanopore: str, r1: str, r2: str, lineage: str, threads: int
):
    drafts_filtered = handle_drafts(drafts)

    ares_out = Path("ares_out")
    try:
        ares_out.mkdir()
    except Exception:
        typer.echo(f"{ares_out} already exists add --f to force overwrite")
        exit()

    pipeline = PolishPipeline(
        root_dir=ares_out.as_posix(),
        drafts=drafts_filtered,
        long_reads=nanopore,
        short_reads=[r1, r2],
        threads=threads,
        lineage=lineage,
    )

    typer.echo("\n****** Running ares ******\n")
    pipeline.run()


@app.command("pilon")
def pilon_command(
    outdir: Path,
    draft: Path,
    r1: str,
    r2: str,
    lineage: str,
    threads: int,
    rounds: int = 4,
):
    if not draft.is_file():
        raise Exception(f"{draft} does not exist")

    try:
        outdir.mkdir()
    except Exception:
        typer.echo(f"{outdir} already exists add --f to force overwrite")
        exit()

    runner = ShortReadPolishRunner(
        root_dir=outdir.as_posix(),
        draft=draft.as_posix(),
        r1=r1,
        r2=r2,
        lineage=lineage,
        threads=threads,
        rounds=rounds,
    )
    best = runner.run(pilon)
    summarize_busco_runs(
        outdir=outdir.as_posix(),
        best=best,
        busco_results=runner.all_busco_runs,
    )


if __name__ == "__main__":
    app()
