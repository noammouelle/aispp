"""Command line: ``python -m detfig configs/generic.yaml [-o out] [--set key=value ...]``."""

import argparse
import sys

from .build import write
from .config import ConfigError, deep_merge, load, parse_override


def main(argv=None):
    ap = argparse.ArgumentParser(prog="detfig", description=__doc__)
    ap.add_argument("configs", nargs="+", help="YAML files, merged left to right on top of the defaults")
    ap.add_argument("-o", "--outdir", default="out", help="output directory (default: out)")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help="override one value, e.g. --set instrument.baseline=10")
    ap.add_argument("--tex-only", action="store_true", help="write .tex but do not run pdflatex")
    args = ap.parse_args(argv)
    try:
        overrides = {}
        for item in args.set:
            overrides = deep_merge(overrides, parse_override(item))
        cfg = load(*args.configs, overrides=overrides)
        files = write(cfg, args.outdir, compile_pdf=not args.tex_only)
    except (ConfigError, ValueError, KeyError, RuntimeError) as err:
        print(f"detfig: error: {err}", file=sys.stderr)
        return 1
    for kind, path in files.items():
        print(f"{kind:>4}: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
