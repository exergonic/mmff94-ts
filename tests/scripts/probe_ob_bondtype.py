"""Probe: OpenBabel's GetBondType rule vs the new reference's bond classes.

For each molecule, read the .mmd block, get OB's per-bond order and
IsAromatic(), apply the OB GetBondType logic with the repo's flag table,
and compare against the FF CLASS column extracted from the opti log.

Run: cd tests/scripts && env -u PYTHONPATH uv run python probe_ob_bondtype.py
"""
import os
import re
import sys
import json

os.environ["BABEL_DATADIR"] = os.path.join(
    os.path.dirname(__import__("openbabel").__file__), "bin", "data")
from openbabel import openbabel as ob

REPO = r"C:/Users/mccan/Code/mmff94-ts"
SUITE_MMD = os.path.join(REPO, "tests", "fixtures", "validation-suite", "MMFF94.mmd")
OPTI_LOG = r"C:/Users/mccan/AppData/Local/Temp/mmff94-761/MMFF94_opti.log"

# our flag table: type -> (arom, sbmb)
FLAGS = {}
props = open(os.path.join(REPO, "src", "mmff94", "parameters",
                          "atom-type-properties.ts"), encoding="utf-8").read()
for m in re.finditer(r"(\d+): \{ crd: \d+, val: \d+, pilp: \d+, mltb: \d+, arom: (\d), lin: \d+, sbmb: (\d)", props):
    FLAGS[int(m.group(1))] = (int(m.group(2)), int(m.group(3)))

def ob_bondtype(mol, i, j):
    """Mirror OB's GetBondType (forcefieldmmff94.cpp, part V p. 620)."""
    bond = mol.GetBond(i, j)
    if bond is None:
        return 0
    if bond.GetBondOrder() != 1 or bond.IsAromatic():
        return 0
    a, b = mol.GetAtom(i + 1), mol.GetAtom(j + 1)
    ta = int(a.GetData("FFAtomType").GetValue())
    tb = int(b.GetData("FFAtomType").GetValue())
    if ta in FLAGS and tb in FLAGS:
        if FLAGS[ta][0] and FLAGS[tb][0]:
            return 1
        if FLAGS[ta][1] and FLAGS[tb][1]:
            return 1
    return 0

# reference classes from the opti log
def ref_bonds(code):
    txt = open(OPTI_LOG, encoding="utf-8", errors="replace").read()
    i = txt.find("New Structure Name/Conformational Index: " + code)
    sec = txt[i:i + 40000]
    j = sec.find("OPTIMOL-ANALYZE>  # bonds")
    k = sec.find("OPTIMOL-ANALYZE>  # ", j + 10)
    out = {}
    for line in sec[j:k].split("\n"):
        m = re.match(r"^ ([A-Za-z0-9]+) #(\d+)\s+([A-Za-z0-9]+) #(\d+)\s+(\d+)\s+(\d+)\s+(\d+)", line)
        if m:
            out[(int(m.group(2)) - 1, int(m.group(4)) - 1)] = int(m.group(7))
    return out

conv = ob.OBConversion()
conv.SetInFormat("mmd")
txt = open(SUITE_MMD, encoding="utf-8").read()
headers = list(re.finditer(r"^\s*\d+\s+\[[A-Za-z0-9]+", txt, re.M))

CODES = ["CEWYIM30", "DAKCEX", "GIGCEE", "KEPKIZ", "SAKGUG", "TAPJUP", "VEWZOM"]
for code in CODES:
    h = next((m for m in headers if "[" + code in m.group(0)), None)
    if not h:
        print(code, "NOT FOUND"); continue
    end = headers[headers.index(h) + 1].start() if headers.index(h) + 1 < len(headers) else len(txt)
    mol = ob.OBMol()
    if not conv.ReadString(mol, txt[h.start():end]):
        print(code, "unreadable"); continue
    ff = ob.OBForceField.FindForceField("MMFF94")
    if not ff.Setup(mol) or not ff.GetAtomTypes(mol):
        print(code, "setup failed"); continue
    ref = ref_bonds(code)
    mism = []
    for (i, j), cls in sorted(ref.items()):
        if i > j: i, j = j, i
        b = mol.GetBond(i + 1, j + 1)
        got = ob_bondtype(mol, i, j)
        # map to the reference class space: double/triple bonds are 0
        if got != cls:
            ti = int(mol.GetAtom(i + 1).GetData("FFAtomType").GetValue())
            tj = int(mol.GetAtom(j + 1).GetData("FFAtomType").GetValue())
            mism.append((i + 1, j + 1, b.GetBondOrder(), b.IsAromatic(),
                         FLAGS.get(ti), FLAGS.get(tj), ti, tj, got, cls))
    print(f"{code}: {len(ref)} bonds, {len(mism)} mismatches vs OB-rule")
    for mm in mism:
        print(f"   bond {mm[0]}-{mm[1]} order {mm[2]} IsArom {mm[3]} "
              f"flags {mm[4]}/{mm[5]} types {mm[6]}/{mm[7]} OB-rule {mm[8]} ref {mm[9]}")
