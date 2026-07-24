from openbabel import openbabel as ob
from openbabel import pybel
import sys
import os

# Locate the openbabel data directory, which contains mmff94.ff
_script_dir = os.path.dirname(os.path.abspath(__file__))

# Try new location first (tests/scripts/.venv/...), then legacy (ob_runs/.venv/...)
_venv_candidates = [
    os.path.join(_script_dir, ".venv"),
    os.path.join(_script_dir, "..", "..", "ob_runs", ".venv"),
]

_ob_data_dir = None
for _candidate in _venv_candidates:
    _test_path = os.path.join(_candidate, "Lib", "site-packages", "openbabel", "bin", "data")
    if os.path.isdir(_test_path):
        _ob_data_dir = _test_path
        break

if _ob_data_dir:
    os.environ['BABEL_DATADIR'] = _ob_data_dir
else:
    print("Warning: Could not find openbabel data directory. Set BABEL_DATADIR manually.")

def get_mmff94_breakdown(input_file, output_file):
    # Read the molecule
    try:
        mol = next(pybel.readfile(input_file.split('.')[-1], input_file))
    except StopIteration:
        print(f"Error: Could not read {input_file}")
        return

    obmol = mol.OBMol
    ff = ob.OBForceField.FindForceField("mmff94")

    if not ff.Setup(obmol):
        print("Error: Could not setup MMFF94 force field.")
        return

    # 1. Extract Energy Components
    e_bond = ff.E_Bond()
    e_angle = ff.E_Angle()
    e_strbnd = ff.E_StrBnd()
    e_torsion = ff.E_Torsion()
    e_oop = ff.E_OOP()
    e_vdw = ff.E_VDW()
    e_elec = ff.E_Electrostatic()
    e_total = ff.Energy()

    # 2. Extract Atom Types and Charges
    # Note: We must explicitly copy FF data to the OBMol to access it easily
    ff.GetAtomTypes(obmol)
    ff.GetPartialCharges(obmol)

    atom_data = []
    for atom in ob.OBMolAtomIter(obmol):
        idx = atom.GetIdx()
        # Atom Type is stored in "FFAtomType" data field after Setup/GetAtomTypes
        atype = atom.GetData("FFAtomType")
        atype_str = atype.GetValue() if atype else "N/A"

        # Partial Charge is stored in "FFPartialCharge"
        pcharge_data = atom.GetData("FFPartialCharge")
        pcharge = float(pcharge_data.GetValue()) if pcharge_data else 0.0

        # Formal Charge
        fcharge_data = atom.GetData("FFFormalCharge")
        fcharge = fcharge_data.GetValue() if fcharge_data else atom.GetFormalCharge()

        # Determine Ring status (Simple check)
        is_ring = "YES" if atom.IsInRing() else "NO"

        atom_data.append({
            "idx": idx,
            "type": atype_str,
            "ring": is_ring,
            "formal_charge": fcharge,
            "partial_charge": pcharge
        })

    # 3. Write to Output File
    with open(output_file, 'w') as f:
        f.write("A T O M   T Y P E S\n\n")
        f.write(f"{'IDX':<8}{'TYPE':<8}{'RING':<8}\n")
        for atom in atom_data:
            f.write(f"{atom['idx']:<8}{atom['type']:<8}{atom['ring']:<8}\n")

        f.write("\nF O R M A L   C H A R G E S\n\n")
        f.write(f"{'IDX':<8}{'CHARGE':<12}\n")
        for atom in atom_data:
            f.write(f"{atom['idx']:<8}{atom['formal_charge']:<12.6f}\n")

        f.write("\nP A R T I A L   C H A R G E S\n\n")
        f.write(f"{'IDX':<8}{'CHARGE':<12}\n")
        for atom in atom_data:
            f.write(f"{atom['idx']:<8}{atom['partial_charge']:<12.6f}\n")

        f.write("\nE N E R G Y\n\n")
        f.write(f"     TOTAL BOND STRETCHING ENERGY =  {e_bond:8.5f} kcal/mol\n")
        f.write(f"     TOTAL ANGLE BENDING ENERGY =  {e_angle:8.5f} kcal/mol\n")
        f.write(f"     TOTAL STRETCH BENDING ENERGY = {e_strbnd:8.5f} kcal/mol\n")
        f.write(f"     TOTAL TORSIONAL ENERGY = {e_torsion:8.5f} kcal/mol\n")
        f.write(f"     TOTAL OUT-OF-PLANE BENDING ENERGY =  {e_oop:8.5f} kcal/mol\n")
        f.write(f"     TOTAL VAN DER WAALS ENERGY =  {e_vdw:8.5f} kcal/mol\n")
        f.write(f"     TOTAL ELECTROSTATIC ENERGY =  {e_elec:8.5f} kcal/mol\n\n")
        f.write(f"TOTAL ENERGY = {e_total:8.5f} kcal/mol\n")

    print(f"Breakdown successfully written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python get_mmff94_breakdown.py <input.sdf> <output.txt>")
    else:
        get_mmff94_breakdown(sys.argv[1], sys.argv[2])
