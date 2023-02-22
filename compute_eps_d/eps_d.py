import numpy as np

seq = {
    'R': 'ARG', 'H': 'HIS', 'K': 'LYS', 'D': 'ASP', 'E': 'GLU',
    'S': 'SER', 'T': 'THR', 'N': 'ASN', 'Q': 'GLN', 'C': 'CYS',
    'U': 'SEC', 'G': 'GLY', 'P': 'PRO', 'A': 'ALA', 'V': 'VAL',
    'I': 'ILE', 'L': 'LEU', 'M': 'MET', 'F': 'PHE', 'Y': 'TYR',
    'W': 'TRP'
}

# Read one letter amino acid sequence from file
filein = 'TDP-43_CTD.dat'
fileout = f"{filein.split('.')[0]}_seq3.dat"
nline = 1
count = 0
with open(fileout, 'w') as fout, open(filein, 'r') as fid:
    for i in fid:
        if not i.startswith('#'):
            for j in i:
                if j in seq:
                    fout.write(f' {seq[j]}')
                    count += 1
                    if count == nline:
                        fout.write('\n')
                        count = 0

# Read sequence and force field parameters
ff_para = 'eps_d_i_i+4.txt'
aalist = {}
with open(ff_para, 'r') as fid:
    for i in fid:
        if not i.startswith('#'):
            tmp = i.split()
            aalist[tmp[0]] = np.loadtxt(tmp[1:], dtype=float)
aakeys = list(aalist.keys())

# This translates each amino acid type into a number
aadih = [aalist[i] for i in aakeys]

# Translate the entire sequence into a number code according to the order in 'aakeys'
chain_id = [aakeys.index(i.split()[0]) for i in open(fileout, 'r')]
chain_length = len(chain_id)
ndih = chain_length - 3
dihedral_pairs = np.array([list(range(i, i + 4)) for i in range(ndih)])
types = list(map(str, range(ndih)))

dih_param = [aalist[i.split()[0]] for i in open(fileout, 'r')]

# Mixing rule
dihed_params = []
for n, i in enumerate(dihedral_pairs):
    if n == 0:
        dihed_params.append([n, f"{(dih_param[i[0]]+dih_param[i[3]]+dih_param[i[3]+1])/3:.6f}"])
    elif n == ndih - 1:
        dihed_params.append([n, f"{(dih_param[i[0]-1]+dih_param[i[0]]+dih_param[i[3]])/3:.6f}"])
    else:
        dihed_params.append([n, f"{(dih_param[i[0]-1]+dih_param[i[0]]+dih_param[i[3]]+dih_param[i[3]+1])/4:.6f}"])

eps_d_all = [item[1] for item in dihed_params]
eps_d = []
n = 1
for i, x in enumerate(eps_d_all):
    if i == eps_d_all.index(x):
        eps_d.append([n,x])
        n=n+1
print(eps_d)
