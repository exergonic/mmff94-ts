"""Generate canonical MMFF94 atom-type references for the validation suite.

Reads the BatchMin-format MMFF94.mmd structure file with OpenBabel and
extracts each molecule's MMFF94 atom types via the force field's
FFAtomType data (OBForceFieldMMFF94::SetTypes assigns the canonical
numeric types; OBForceField::GetAtomTypes copies them onto the input
molecule as per-atom "FFAtomType" data).

Output: JSON keyed by molecule code -> list of atom type integers.
Molecules whose force-field setup fails (metals, untyped elements) are
recorded in a "skipped" list and excluded from the reference.

Usage: uv run --project tests/scripts python tests/scripts/extract_suite_types.py
"""

import json
import os
import re
import sys

from openbabel import openbabel as ob

# openbabel-wheel on Windows points BABEL_DATADIR at share/openbabel/3.1.0
# (which only holds splash.png); the parameter files live in bin/data.
os.environ["BABEL_DATADIR"] = os.path.join(
    os.path.dirname(ob.__file__), "bin", "data")

SUITE_MMD = os.path.join(os.path.dirname(__file__), "..", "fixtures",
                         "validation-suite", "MMFF94.mmd")
OUT_JSON = os.path.join(os.path.dirname(__file__), "..", "fixtures",
                        "validation-suite", "mmff94-atom-types.json")

# Each molecule in the .mmd file starts with a header line of the form
# "   <atom count>   [CODE,...]". Split on those headers so each molecule
# can be fed to OpenBabel individually: the binding's ReadString is
# one-shot (it always re-reads the first molecule from the string), so a
# whole-file read cannot be iterated.
HEADER_RE = re.compile(r"^\s*\d+\s+\[[A-Za-z0-9]+", re.M)

conv = ob.OBConversion()
if not conv.SetInFormat("mmd"):
    print("error: openbabel has no mmd reader", file=sys.stderr)
    sys.exit(1)

ff = ob.OBForceField.FindForceField("MMFF94")

reference: dict[str, list[int]] = {}
skipped: dict[str, str] = {}

with open(SUITE_MMD) as f:
    content = f.read()

headers = list(HEADER_RE.finditer(content))
mol = ob.OBMol()
for i, header in enumerate(headers):
    end = headers[i + 1].start() if i + 1 < len(headers) else len(content)
    block = content[header.start():end]

    if not conv.ReadString(mol, block):
        skipped[f"mol-{i}"] = "unreadable block"
        continue

    # OpenBabel titles come back as "[CODE,10,10,S,k]" — strip to the code.
    code = mol.GetTitle().split(",")[0].strip(" []")

    if not ff.Setup(mol):
        skipped[code] = "force-field setup failed"
        continue
    if not ff.GetAtomTypes(mol):
        skipped[code] = "GetAtomTypes failed"
        continue
    types = []
    ok = True
    for a in range(mol.NumAtoms()):
        dp = mol.GetAtom(a + 1).GetData("FFAtomType")
        if dp is None:
            ok = False
            break
        types.append(int(dp.GetValue()))
    if not ok:
        skipped[code] = "missing FFAtomType data"
        continue
    reference[code] = types

with open(OUT_JSON, "w") as f:
    json.dump({"molecules": reference, "skipped": skipped}, f, indent=1)

print(f"{len(reference)} molecules typed, {len(skipped)} skipped, "
      f"{len(headers)} blocks found")
for code, why in sorted(skipped.items()):
    print(f"  skip {code}: {why}")
