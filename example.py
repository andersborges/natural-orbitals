import numpy as np
import pickle
from NO import NaturalOrbitals as NO
from NO import plot_basis


#### import Hamiltonian/overlap matrix #####
#calculated with gpaw. Dumped with dump_hamiltonian_parallel() as pickles.  
h = np.loadtxt('h.txt') 
s = np.loadtxt('s.txt')

#### import Atoms object #####
from ase.io import read
mol = read('CH4.txt')
ne = 4 # number of electrons pairs (can be found in gpaw output file)

#### generate library which links atoms to basis functions ####
from gpaw.lcao.tools import get_bfi2
basis = 'sz'
basis_dic = {}
for atom in mol:
	basis_dic[atom.index] = get_bfi2(mol.get_chemical_symbols(), basis, [atom.index])


#### Initialize ####
CH4 = NO(h, s,ne, mol, basis_dic, model = True) #initialize

print "Natural electronic contribution to atomic charges: ", CH4.get_natural_charges() # natural atomic charges (nuclear charge not included)

#get natural hybrids
NHOs = CH4.get_natural_hybrids()
#plot natural hybrids using gpaw
plot_basis(NHOs, mol, range(8), basis = 'sz', folder_name='./NHOs',)

#plot pre-natural hybrids using gpaw
plot_basis(CH4.PNHO, mol, range(8), basis = 'sz', folder_name='./PNHOs',)


