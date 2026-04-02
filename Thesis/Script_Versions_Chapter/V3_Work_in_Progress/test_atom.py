#!/usr/bin/python3


import pandas as pd
import io
import itertools

import Atom


################################################################################
################################################################################
################################################################################
#header = ['opls_name', 'name', 'atom_num', 'mass', 'charge', 'ptype', 'sigma', 'epsilon']
#header = ['opls_name', 'name', 'at.num', 'mass', 'charge', 'ptype', 'sigma', 'epsilon']
header = ['type', 'name', 'at.num', 'mass', 'charge', 'ptype', 'sigma', 'epsilon']
df_opls = pd.read_csv('ffnonbonded.itp', skiprows=3, sep=r'\s+', comment=';',
                      header=None, names=header,
                      index_col='type')
#df.set_index('opls_name')
print(df_opls.head())
print(df_opls.loc['opls_002'])
print(df_opls.loc['opls_002']['epsilon'])
#exit(0)

#for index, row in df.iterrows():
#    print(row['opls_name'])


# read topology file *.itp

f = open('EthaneDiol.itp')
header_to_end = itertools.dropwhile(lambda x: x.strip() != '[ atoms ]', f)

atoms_section = itertools.takewhile(lambda x: x.strip() != '', header_to_end)
data_str = ''.join(atoms_section)

#for line in atoms_section:
#    print(line)

f.close()

print(data_str)
#header = ['idx', 'opls_name', 'num1', 'res_name', 'two_letter', 'res_num', 'charge', 'mass']
header = ['nr', 'type', 'resnr', 'residue', 'atom', 'cgnr', 'charge', 'mass']
df = pd.read_csv(io.StringIO(data_str), skiprows=[0],
                 sep=r'\s+', header=None, names=header, comment=';')
#df.columns = ['nr', 'type', 'resnr', 'residue', 'atom', 'cgnr', 'charge', 'mass']
print(df)    

atom_list = []
for index, row in df.iterrows():
    atom = Atom.Atom()
    type_str = row['type']
    print(type_str)
    if (type_str in df_opls.index):
        row_opls = df_opls.loc[type_str]
        for key, value in row_opls.items():
            atom[key] = value
    for key, value in row.items():
        print(key, value)
        atom[key] = value
    atom_list.append(atom)

for atom in atom_list:
    print(atom['nr'], atom['type'], atom['ptype'], atom['charge'], atom['sigma'], atom['epsilon'])
    
exit(0)


#for index, row in df.iterrows():
#    atom = Atom.Atom()
#    for key, value in row.items():
#        atom[key] = value
#    atom_list.append(atom)
print(df)



################################################################################
################################################################################
################################################################################
