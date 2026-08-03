"""Probe OpenBabel MMFF94 typing on synthetic SMILES molecules to confirm
the phosphate/sulfonate oxygen rules generalize beyond the validation suite."""
import os
import sys

from openbabel import openbabel as ob

os.environ["BABEL_DATADIR"] = os.path.join(
    os.path.dirname(ob.__file__), "bin", "data")

SMILES = [
    ("sulfite", "O=S([O-])[O-]"),
    ("bisulfite", "O=S(O)[O-]"),
    ("sulfate", "O=S(=O)([O-])[O-]"),
    ("bisulfate", "O=S(=O)(O)[O-]"),
    ("methanesulfonate", "CS(=O)(=O)[O-]"),
    ("methanesulfonic acid", "CS(=O)(=O)O"),
    ("dimethyl sulfone", "CS(C)(=O)=O"),
    ("dimethyl sulfoxide", "CS(C)=O"),
    ("methanesulfinic acid", "CS(=O)O"),
    ("methanesulfinamide", "CS(=O)N"),
    ("sulfamide", "NS(=O)(=O)N"),
    ("methanesulfonamide", "CS(=O)(=O)N"),
    ("dimethyl sulfate", "CO[S](=O)(=O)OC"),
    ("methyl hydrogen sulfate", "COS(=O)(=O)O"),
    ("thiosulfate", "O=S(=O)([O-])[S-]"),
    ("sulfine", "C=S=O"),
    ("sulfoximine", "CS(=O)(=N)C"),
    ("sulfonimidamide", "CS(=O)(=N)N"),
    ("methanesulfinate anion", "CS(=O)[O-]"),
    ("sulfite monoester anion", "CO[S](=O)[O-]"),
    ("sulfurous acid", "OS(O)=O"),
    ("hypochlorous acid", "OCl"),
    ("hypochlorite", "[O-]Cl"),
    ("chlorate", "Cl(=O)(=O)[O-]"),
    ("thiocarbonyl", "CC(=S)C"),
    ("thioacetamide", "CC(=S)N"),
    ("thiourea", "NC(=S)N"),
    ("dithioacetate", "CC(=S)[S-]"),
    ("thioformaldehyde", "C=S"),
    ("phosphoric acid", "OP(=O)(O)O"),
    ("dihydrogen phosphate", "OP(=O)(O)[O-]"),
    ("hydrogen phosphate", "OP(=O)([O-])[O-]"),
    ("phosphate", "O=P([O-])([O-])[O-]"),
    ("trimethyl phosphate", "COP(=O)(OC)OC"),
    ("dimethyl phosphate", "COP(=O)(OC)[O-]"),
    ("phosphine oxide", "OP"),
    ("trimethylphosphine oxide", "CP(C)(C)=O"),
    ("phosphine", "P"),
    ("trimethylphosphine", "CP(C)C"),
    ("phosphorous acid", "OP(O)O"),
    ("phosphoramidic acid", "NP(=O)(O)O"),
    ("phosphine sulfide", "CP(C)(C)=S"),
    ("methylphosphonic acid", "CP(=O)(O)O"),
    ("perchlorate", "Cl(=O)(=O)(=O)[O-]"),
    ("hypochlorite", "[O-]Cl"),
    ("chlorate", "Cl(=O)(=O)[O-]"),
]

conv = ob.OBConversion()
conv.SetInFormat("smi")
ff = ob.OBForceField.FindForceField("MMFF94")
mol = ob.OBMol()

for name, smi in SMILES:
    mol.Clear()
    if not conv.ReadString(mol, smi):
        print(f"{name:28s} {smi:28s} UNREADABLE")
        continue
    mol.AddHydrogens()
    if not ff.Setup(mol):
        print(f"{name:28s} {smi:28s} SETUP FAILED")
        continue
    ff.GetAtomTypes(mol)
    # map: atom -> (element, type)
    parts = []
    for a in range(mol.NumAtoms()):
        atom = mol.GetAtom(a + 1)
        dp = atom.GetData("FFAtomType")
        t = int(dp.GetValue()) if dp else -1
        parts.append(f"{atom.GetAtomicNum()}:{t}")
    print(f"{name:28s} {smi:28s} {' '.join(parts)}")
