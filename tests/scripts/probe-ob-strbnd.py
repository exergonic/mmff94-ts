#!/usr/bin/env python
"""Capture OpenBabel's high-verbosity energy log for a suite molecule."""
import os
import sys
import tempfile

from openbabel import openbabel as ob

os.environ['BABEL_DATADIR'] = os.path.join(os.path.dirname(ob.__file__), 'bin', 'data')
name = sys.argv[1]
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

fd, path = tempfile.mkstemp()
saved_out = os.dup(1)
saved_err = os.dup(2)
os.dup2(fd, 1)
os.dup2(fd, 2)
ob.obErrorLog.SetOutputLevel(ob.OBFF_LOGLVL_HIGH)
try:
    ob.obErrorLog.SetOutputStream(sys.stdout)
except TypeError:
    pass
e = ff.Energy()
os.dup2(saved_out, 1)
os.dup2(saved_err, 2)
os.close(fd)
data = open(path, encoding='utf-8', errors='replace').read()
os.unlink(path)
print(f'energy={e} bytes={len(data)}')
for line in data.splitlines():
    print(line)
