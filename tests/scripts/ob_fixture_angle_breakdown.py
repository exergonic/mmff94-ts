#!/usr/bin/env python
"""OpenBabel per-angle MMFF94 breakdown for an SDF fixture.

Reads a fixture SDF, assigns MMFF94 types + partial charges, and prints
the HIGH-verbosity ANGLE BENDING interaction log (the per-angle
energies with the atoms involved), plus the per-term totals.

Usage (from tests/scripts, with the uv project's venv):
  uv run ob_fixture_angle_breakdown.py NAME   (NAME = SDF stem, e.g. trimethylamine)
"""
import os
import re
import sys
import tempfile

from openbabel import openbabel as ob

os.environ["BABEL_DATADIR"] = os.path.join(os.path.dirname(ob.__file__), "bin", "data")

SDF_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fixtures", "sdf")


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "trimethylamine"
    conv = ob.OBConversion()
    conv.SetInFormat("sdf")
    mol = ob.OBMol()
    ok = conv.ReadFile(mol, f"{SDF_DIR}/{name}.sdf")
    if not ok:
        print(f"cannot read {name}.sdf")
        sys.exit(1)

    ff = ob.OBForceField.FindForceField("MMFF94")
    ff.Setup(mol)
    ff.GetAtomTypes(mol)
    ff.GetPartialCharges(mol)
    types = [mol.GetAtom(i + 1).GetType() for i in range(mol.NumAtoms())]

    fd, path = tempfile.mkstemp()
    saved_err = os.dup(2)
    os.dup2(fd, 2)
    ff.SetLogToStdErr()
    ff.SetLogLevel(ob.OBFF_LOGLVL_HIGH)
    total = ff.Energy()
    os.dup2(saved_err, 2)
    os.close(fd)
    log = open(path, encoding="utf-8", errors="replace").read()
    os.unlink(path)

    print(f"{name}: {mol.NumAtoms()} atoms")
    print("  types:   ", " ".join(types))
    print(f"  total: {total:.5f}")

    lines = log.splitlines()
    # The angle section: from the de-spaced header to its TOTAL line.
    start = None
    for i, line in enumerate(lines):
        if line.strip() and line.strip()[0].isalpha() and "ANGLEBENDING" in line.replace(" ", "").upper():
            start = i
            break
    if start is None:
        print("  (no angle section in log)")
        return
    for l in lines[start:]:
        print(l)
        if "TOTAL" in l and "ENERGY =" in l:
            break


if __name__ == "__main__":
    main()
