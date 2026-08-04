#!/usr/bin/env python
"""OpenBabel MMFF94 workflow exemplar.

The complete reference-computation loop used across this repo's
analysis scripts, in one readable pass:

  1. the BABEL_DATADIR quirk (openbabel-wheel on Windows)
  2. reading a molecule (a validation-suite mmd block by default)
  3. MMFF94 setup with canonical atom types and partial charges
  4. the energy with a per-component breakdown
  5. the HIGH-verbosity per-interaction log (the debugging workhorse —
     `--verbose strbnd` prints every stretch-bend interaction)

Usage (run from the repo root, with tests/scripts/.venv):
  tests/scripts/.venv/Scripts/python tests/scripts/ob_energy_breakdown.py [NAME] [--verbose [TERM]]
  e.g. COYVIV (delocalized-N sulfone), FUCTIG01 (nitrate), VIYPAU (acetal)

The force field itself is not modified here; the only environment
adjustment is BABEL_DATADIR, which every script in this directory
needs on Windows.
"""

import os
import re
import sys
import tempfile

from openbabel import openbabel as ob

# openbabel-wheel on Windows points BABEL_DATADIR at share/openbabel/3.1.0,
# which holds only splash.png. The parameter files (mmff*.par, types.txt)
# live in bin/data, and without them MMFF94's Setup() fails — so every
# script here redirects the data dir to the wheel's actual data folder.
os.environ["BABEL_DATADIR"] = os.path.join(os.path.dirname(ob.__file__), "bin", "data")

SUITE_MMD = "tests/fixtures/validation-suite/MMFF94.mmd"


def load_suite_molecule(code):
    """Read one molecule's block from the suite mmd and parse it with OB.

    The mmd file is a concatenation of per-molecule blocks; each starts
    with a header line "<natoms> [CODE,...]". OpenBabel reads a block as
    its own OBMol through the mmd format reader.
    """
    mmd = open(SUITE_MMD, encoding="utf-8", errors="replace").read()
    blocks, cur = [], None
    for line in mmd.splitlines():
        m = re.match(r"^\s*\d+\s+\[(\w+),", line)
        if m:
            if cur is not None:
                blocks.append(cur)
            cur = [line]
        elif cur is not None:
            cur.append(line)
    if cur is not None:
        blocks.append(cur)

    block = next((b for b in blocks if f"[{code}," in b[0]), None)
    if block is None:
        sys.exit(f"molecule {code} not found in {SUITE_MMD}")

    conv = ob.OBConversion()
    conv.SetInFormat("mmd")
    mol = ob.OBMol()
    if not conv.ReadString(mol, "\n".join(block) + "\n"):
        sys.exit(f"OpenBabel could not read {code}")
    return mol


def setup_mmff94(mol):
    """Find the MMFF94 force field and set it up on the molecule.

    Setup() assigns the MMFF94 atom types and collects every interaction
    (bonds, angles, torsions, vdW pairs...). Failure means untyped atoms
    or missing parameters — the error log names them.
    """
    ff = ob.OBForceField.FindForceField("MMFF94")
    if not ff.Setup(mol):
        sys.exit("MMFF94 setup failed (see the warnings above)")
    return ff


def compute_energy(ff, verbose=False, term=""):
    """Compute the energy, capturing the per-component and per-interaction log.

    The Python binding has no GetLog(): the force field writes its
    verbose log through the C++ streams. We redirect fd 2 to a temp
    file and ask for stderr output — std::cerr is unbuffered, so every
    write lands in the file immediately. (stdout would NOT work:
    std::cout is buffered, the file comes out empty and the whole log
    leaks to the terminal when the buffer flushes at exit.)
    """
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

    # Per-component totals: "TOTAL BOND STRETCHING ENERGY =  0.16767 kcal/mol"
    # (double space after the equals sign; hyphenated names like
    # OUT-OF-PLANE).
    components = {}
    for line in log.splitlines():
        m = re.match(r"\s*TOTAL ([A-Z][A-Z -]*) ENERGY =\s+([-\d.]+)", line)
        if m:
            components[m.group(1).strip()] = float(m.group(2))

    if verbose:
        for line in log.splitlines():
            if not term or term.lower() in line.lower():
                print(line)
    return total, components


def main():
    code = sys.argv[1] if len(sys.argv) > 1 else "COYVIV"
    verbose = "--verbose" in sys.argv[1:]
    term = next((a for a in sys.argv[2:] if not a.startswith("--") and a != code), "")

    mol = load_suite_molecule(code)
    ff = setup_mmff94(mol)

    # Canonical atom types and partial charges, copied onto the molecule
    # (both methods take the molecule: they fill its per-atom data).
    ff.GetAtomTypes(mol)
    types = [mol.GetAtom(i + 1).GetType() for i in range(mol.NumAtoms())]
    ff.GetPartialCharges(mol)
    charges = [mol.GetAtom(i + 1).GetPartialCharge() for i in range(mol.NumAtoms())]

    total, components = compute_energy(ff, verbose, term)

    print(f"{code}: {mol.NumAtoms()} atoms")
    print("  types:   ", " ".join(types))
    print("  charges: ", " ".join(f"{c:+.3f}" for c in charges))
    print(f"  total energy: {total:.4f} kcal/mol")
    for name, val in components.items():
        print(f"  {name:24s} {val:10.4f}")


if __name__ == "__main__":
    main()
