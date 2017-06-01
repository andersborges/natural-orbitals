import numpy as np
import pickle
from NO import NaturalOrbitals as NO
from NO import plot_basis


# Import Hamiltonian/overlap matrix. This example has been 
# calculated/optimized with GPAW. H and S can be saved to disk 
# with dump_hamiltonian_parallel() as pickles and be
# converted to .txt-files.   
h = np.loadtxt('files_CH4/h.txt') 
s = np.loadtxt('files_CH4/s.txt')

# Import Atoms object. 
from ase.io import read
mol = read('files_CH4/CH4.txt')
ne = 4 # Set number of electrons pairs (can be found in gpaw output file or from the calculator object). 

# Generate library which links atoms to basis functions (can also be done using get_bfi2 from gpaw)
basis_dic = {0:range(4), 1:[4], 2:[5], 3:[6], 4:[7]}

# Initialize 
CH4 = NO(h=h, s=s,ne=ne, mol=mol, basis_dic=basis_dic, model = './files_CH4/model') # model requires ase installed. 

print "Natural electronic contribution to atomic charges: ", CH4.get_natural_charges() # natural atomic charges (nuclear charge not included)
#get natural hybrids
NHOs = CH4.get_natural_bonding_orbitals() # generates 

#Plot Natural Hybrid Orbital (NHO) using GPAW. Requires gpaw sz basis set is available.
# Basis can be generated with "gpaw-basis -t sz H C". Cube files can be found in files_CH4/NHOs.
plot_basis(NHOs, mol, range(8), basis = 'sz', folder_name='./files_CH4/NHOs')

#plot pre-natural hybrids using gpaw
plot_basis(CH4.PNHO, mol, range(8), basis = 'sz', folder_name='./files_CH4/PNHOs')


