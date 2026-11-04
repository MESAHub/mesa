"""Render both model variants and refresh the four cases' README figures."""

from pathlib import Path
import subprocess
import tempfile
import shutil

source = Path(__file__).resolve().with_name("main.tex")
models = source.parent.parent
mesa = models.parents[2]
names = ("run_star_extras", "run_star_extras_layout")

for case, binary in (
    ("dev_TDC_single_star_ppisn", False),
    ("dev_TDC_ppisn_from_binary", True),
):
    docs = models / case / "docs"
    docs.mkdir(exist_ok=True)
    definition = r"\def\BinaryProduct{1}" if binary else ""
    tex = definition + r"\input{" + source.as_posix() + "}"
    with tempfile.TemporaryDirectory() as work:
        subprocess.run(
            [
                "pdflatex",
                "-halt-on-error",
                "-interaction=nonstopmode",
                "-jobname=main",
                "-output-directory",
                work,
                tex,
            ],
            cwd=work,
            check=True,
        )
        for page, name in enumerate(names, start=1):
            subprocess.run(
                [
                    "pdftoppm",
                    "-f",
                    str(page),
                    "-l",
                    str(page),
                    "-singlefile",
                    "-r",
                    "180",
                    "-png",
                    str(Path(work) / "main.pdf"),
                    str(docs / name),
                ],
                check=True,
            )
            print(docs / (name + ".png"))

for case in ("ppisn", "pisn"):
    docs = mesa / "star/test_suite" / case / "docs"
    docs.mkdir(exist_ok=True)
    for name in names:
        shutil.copyfile(
            models / "dev_TDC_single_star_ppisn/docs" / (name + ".png"),
            docs / (name + ".png"),
        )
        print(docs / (name + ".png"))
