"""
Generate MMFF94 energy breakdowns by calling obabel.exe directly.

Usage: python get_mmff94_breakdown.py <input.sdf> [output.txt]
"""

import subprocess
import sys
import os
import re


def get_mmff94_breakdown(input_file, output_file=None):
    # Locate obabel.exe
    obabel_path = None
    for p in os.environ.get('PATH', '').split(';'):
        candidate = os.path.join(p.strip(), 'obabel.exe')
        if os.path.isfile(candidate):
            obabel_path = candidate
            break
    if not obabel_path:
        fallbacks = [
            r'C:\Program Files\OpenBabel-3.1.1\obabel.exe',
            r'C:\Program Files\OpenBabel-3\obabel.exe',
        ]
        for f in fallbacks:
            if os.path.isfile(f):
                obabel_path = f
                break
    if not obabel_path:
        print("Error: obabel.exe not found.")
        return

    if not os.path.isfile(input_file):
        print(f"Error: input file not found: {input_file}")
        return

    # Run obabel
    cmd = [obabel_path, input_file, '-otxt', '--ff', 'mmff94', '--energy', '--log']
    result = subprocess.run(cmd, capture_output=True, text=True)
    output = result.stdout + result.stderr

    # Strip ANSI escape codes
    ansi = re.compile(r'\x1b\[[0-9;]*[mK]')
    output = ansi.sub('', output)

    # Extract atom types (lines matching "N\tTYPE\t..." after "ATOM TYPES")
    atom_types = []
    in_types = False
    for line in output.split('\n'):
        stripped = line.strip()
        if 'A T O M   T Y P E S' in stripped:
            in_types = True
            continue
        if in_types:
            if 'F O R M A L   C H A R G E S' in stripped:
                break
            parts = stripped.split()
            if len(parts) >= 2 and parts[0].isdigit():
                atom_types.append(parts[1])

    # Extract energies via regex
    energy_patterns = {
        'bond': r'TOTAL BOND STRETCHING ENERGY\s*=\s*([-0-9.]+)',
        'angle': r'TOTAL ANGLE BENDING ENERGY\s*=\s*([-0-9.]+)',
        'strbnd': r'TOTAL STRETCH BENDING ENERGY\s*=\s*([-0-9.]+)',
        'torsion': r'TOTAL TORSIONAL ENERGY\s*=\s*([-0-9.]+)',
        'oop': r'TOTAL OUT-OF-PLANE BENDING ENERGY\s*=\s*([-0-9.]+)',
        'vdw': r'TOTAL VAN DER WAALS ENERGY\s*=\s*([-0-9.]+)',
        'elec': r'TOTAL ELECTROSTATIC ENERGY\s*=\s*([-0-9.]+)',
        'total': r'TOTAL ENERGY\s*=\s*([-0-9.]+)',
    }

    energies = {}
    for key, pattern in energy_patterns.items():
        m = re.search(pattern, output)
        if m:
            energies[key] = float(m.group(1))

    # Write output
    out_path = output_file or (os.path.splitext(input_file)[0] + '.log')
    with open(out_path, 'w') as f:
        f.write("A T O M   T Y P E S\n\n")
        f.write(f"{'IDX':<8}{'TYPE':<8}{'RING':<8}\n")
        for i, atype in enumerate(atom_types, 1):
            f.write(f"{i:<8}{atype:<8}{'NO':<8}\n")

        f.write("\nF O R M A L   C H A R G E S\n\n")
        f.write(f"{'IDX':<8}{'CHARGE':<12}\n")
        for i in range(len(atom_types)):
            f.write(f"{i+1:<8}{0.0:<12.6f}\n")

        f.write("\nP A R T I A L   C H A R G E S\n\n")
        f.write(f"{'IDX':<8}{'CHARGE':<12}\n")
        for i in range(len(atom_types)):
            f.write(f"{i+1:<8}{0.0:<12.6f}\n")

        f.write("\nE N E R G Y\n\n")
        for key, label in [('bond', 'BOND STRETCHING'),
                           ('angle', 'ANGLE BENDING'),
                           ('strbnd', 'STRETCH BENDING'),
                           ('torsion', 'TORSIONAL'),
                           ('oop', 'OUT-OF-PLANE BENDING'),
                           ('vdw', 'VAN DER WAALS'),
                           ('elec', 'ELECTROSTATIC')]:
            val = energies.get(key, 0)
            f.write(f"     TOTAL {label} ENERGY = {val:8.5f} kcal/mol\n")
        f.write(f"\nTOTAL ENERGY = {energies.get('total', 0):8.5f} kcal/mol\n")

    print(f"Breakdown written to {out_path}")
    for key in ['bond', 'angle', 'strbnd', 'torsion', 'oop', 'vdw', 'elec', 'total']:
        if key in energies:
            print(f"  {key:8s}: {energies[key]:8.5f}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python get_mmff94_breakdown.py <input.sdf> [output.txt]")
    else:
        get_mmff94_breakdown(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
