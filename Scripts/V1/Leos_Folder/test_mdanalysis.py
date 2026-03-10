#!/usr/bin/python3

import MDAnalysis

data = MDAnalysis.Universe('./Ethanol_UA/nvt_vacuum2.gro')
print(data.filename)
print(data.dimensions)
print(data.atoms.n_atoms)
for atom in data.atoms:
    print(atom.name, atom.type, atom.resname, atom.resnum, atom.position)
print(data.atoms)
#print(data.trajectory.get('title'))
