#!/usr/bin/env python
"""Dump OpenBabel's per-interaction energy log for a suite molecule (all terms)."""
import os
import sys
import tempfile

from openbabel import openbabel as ob

os.environ['BABEL_DATADIR'] = os.path.join(os.path.dirname(ob.__file__), 'bin', 'data')
name = sys.argv[1]
term = sys.argv[2] if len(sys.argv) > 2 else ''
mmd = open('tests/fixtures/validation-suite/MMFF94.mmd', encoding='utf-8', errors='replace').read()
blocks, cur = [], None
for line in mmd.splitlines():
    if line.strip() and line.strip()[0].isdigit() and '[' in line:
        if cur is not None:
            blocks.append(cur)
        cur = [line]
    elif cur is not None:
        cur.append(line)
if cur is not None:
    blocks.append(cur)
block = next(b for b in blocks if f'[{name},' in b[0])

conv = ob.OBConversion()
conv.SetInFormat('mmd')
mol = ob.OBMol()
conv.ReadString(mol, '\n'.join(block) + '\n')
ff = ob.OBForceField.FindForceField('MMFF94')
ff.Setup(mol)

fd, path = tempfile.mkstemp(suffix='.log')
saved_out = os.dup(1)
os.dup2(fd, 1)
ff.SetLogToStdOut()
ff.SetLogLevel(ob.OBFF_LOGLVL_HIGH)
e = ff.Energy()
os.dup2(saved_out, 1)
os.close(fd)
data = open(path, encoding='utf-8', errors='replace').read()
os.unlink(path)
print(f'energy={e:.6f}')
for line in data.splitlines():
    low = line.lower()
    if not term or term.lower() in low:
        print(line)
