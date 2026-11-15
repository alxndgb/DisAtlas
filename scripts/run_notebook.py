"""Execute every code cell of a notebook in order (outputs are not written back)."""
import sys, os, nbformat
nb = nbformat.read(sys.argv[1], as_version=4)
os.chdir(os.path.dirname(os.path.abspath(sys.argv[1])))
try:  # notebook cells call the IPython builtin display(); fall back to print outside IPython
    from IPython.display import display
except ImportError:
    display = print
g = {"__name__": "__main__", "display": display}
for i, c in enumerate(nb.cells):
    if c.cell_type == "code":
        src = "\n".join(l for l in c.source.splitlines() if not l.lstrip().startswith(("%", "!")))
        exec(compile(src, f"{sys.argv[1]}[cell {i}]", "exec"), g)
print(f"ok {sys.argv[1]}")
