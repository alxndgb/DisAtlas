"""Execute the exploratory notebooks top to bottom and store their (small, text-only) outputs.

Cells run in-process in one namespace per notebook, so no Jupyter kernel or network port is
needed; every `assert` in a notebook checks a recomputed value against the manuscript. Usage,
from the repository root:

    python notebooks/execute_notebooks.py            # all notebooks
    python notebooks/execute_notebooks.py 04_headline_evaluation.ipynb
"""
import contextlib
import io
import os
import sys
import traceback
import warnings
from pathlib import Path

import nbformat

HERE = Path(__file__).resolve().parent


def run(path):
    nb = nbformat.read(path, as_version=4)
    ns = {"__name__": "__main__"}
    k = 0
    for cell in nb.cells:
        if cell.cell_type != "code":
            continue
        k += 1
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                exec(compile(cell.source, f"{path.name}:cell{k}", "exec"), ns)
        except Exception:
            cell.execution_count, cell.outputs = k, []
            print(buf.getvalue())
            traceback.print_exc()
            return False
        cell.execution_count = k
        text = buf.getvalue()
        cell.outputs = [nbformat.v4.new_output("stream", name="stdout", text=text)] if text else []
    nbformat.write(nb, path)
    return True


def main():
    warnings.filterwarnings("ignore")
    names = sys.argv[1:] or sorted(p.name for p in HERE.glob("[0-9][0-9]_*.ipynb"))
    os.chdir(HERE)
    failed = [n for n in names if not run(HERE / n)]
    for n in names:
        print(f"{n}: {'FAILED' if n in failed else 'ok'}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
