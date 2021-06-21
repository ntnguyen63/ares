import typer

app = typer.Typer()


@app.command()
def ares(drafs: str, nanopore: str, r1: str, r2: str):
    typer.echo(f"{nanopore}")


if __name__ == "__main__":
    app()
