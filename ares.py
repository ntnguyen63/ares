from pathlib import Path

import typer

from ares.polisher import PolishPipeline

app = typer.Typer()


@app.command()
def ares(
    drafts: str, nanopore: str, r1: str, r2: str, busco_lineage: str, t: int
):
    drafts_dir = Path(drafts)
    if not drafts_dir.is_dir():
        raise Exception(f"{drafts} is not a directory")

    drafts_filtered = [
        draft for draft in drafts_dir.glob("*.fasta") if draft.is_file()
    ]
    if not drafts_filtered:
        raise Exception(f"No fasta files in {drafts}")

    ares_out = Path("ares")
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
        threads=t,
        lineage=busco_lineage,
    )

    typer.echo("\n****** Running ares ******\n")
    pipeline.run()


if __name__ == "__main__":
    app()
