import os
from openbabel import openbabel as ob
import openbabel as ob2

os.environ["BABEL_DATADIR"] = os.path.join(os.path.dirname(ob2.__file__), "bin", "data")

conv = ob.OBConversion()
conv.SetInFormat("sdf")
mol = ob.OBMol()
conv.ReadFile(mol, "cswat.sdf")
ff = ob.OBForceField.FindForceField("MMFF94")
ok = ff.Setup(mol)
print("setup:", ok)
print("types:", [a.GetType() for a in ob.OBMolAtomIter(mol)])
print("E total:", ff.Energy())
for name in ["Bond", "Angle", "StretchBend", "Torsion", "OutOfPlane", "VDW", "Electrostatic"]:
    try:
        print(name, getattr(ff, name + "Energy")())
    except Exception as ex:
        print(name, "n/a:", ex)
