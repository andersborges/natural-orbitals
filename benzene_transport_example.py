import numpy as np
from NO import NaturalOrbitals as NO
from NO import plot_basis
from ase.io import read, write
from gpaw.lcao.tools import get_bfi2

# Import Hamiltonian and overlap matrix. This example is calculated using GPAW. 
h = np.loadtxt('files_benzene/h.txt')
s = np.loadtxt('files_benzene/s.txt')

# Import Atoms object 
mol = read('files_benzene/benzene.txt')

# Specify number of electron pairs. 
ne = 15

# Generate library which links atoms to basis functions 
basis_dic = {}
for atom in mol:
	basis_dic[atom.index]=get_bfi2(mol.get_chemical_symbols(), 'dzp', [atom.index])

# Initialize class
benz = NO(h=h, s=s, ne=ne, mol=mol, basis_dic = basis_dic, NHO_model=False, NBO_model=False)

np.set_printoptions(suppress=True, precision =1)
# Get Natural Atomic Orbitals (NAOs)
NAO = benz.get_NAO()

# Get indices of NAOs with occupancy of 1. For benzene, these correspond to the pz-orbitals.
pzs = benz.get_single_occupancy_NAO_indices(threshold=0.01)

# Uncomment the line below to plot pz orbitals. This requires that GPAW is installed. 
#plot_basis(NAO, mol, pzs, basis = 'dzp', folder_name='./files_benzene/pzs') 

##### transport #######
# Specify energy grid. 
energies = np.arange(benz.E_F-4.0,benz.E_F+4.0, 0.01) # Attribute "E_F" is halfway between HOMO and LUMO
# Specify wide-band self energies. Below we assume only one pz couples to each lead. These sites sit "para-" to each other"
SigmaL = np.zeros(h.shape, dtype = 'complex')
SigmaL[pzs[0], pzs[0]] = 0.1j
SigmaR = np.zeros(h.shape, dtype = 'complex')
SigmaR[pzs[3], pzs[3]] = 0.1j

# calculate transmission for para
t_para = benz.get_transmission(energies, SigmaL, SigmaR )

# repeat for meta and ortho
SigmaR = np.zeros(h.shape, dtype = 'complex')
SigmaR[pzs[2], pzs[2]] = 0.1j
t_meta = benz.get_transmission(energies, SigmaL, SigmaR)

SigmaR = np.zeros(h.shape, dtype = 'complex')
SigmaR[pzs[1], pzs[1]] = 0.1j
t_ortho = benz.get_transmission(energies, SigmaL, SigmaR)

#### plot transmission through pi system. ####
from matplotlib import pylab as plt
plt.figure(1, (4,4))
plt.semilogy(energies, t_para, label= 'para')
plt.semilogy(energies, t_meta, label= 'meta')
plt.semilogy(energies, t_ortho, label= 'ortho')

plt.ylim([10**(-5), 2])
plt.xlim([-4, 4])
plt.legend(fontsize=8, loc = 'lower left')
plt.xlabel(r'$E-E_F$ (eV)')
plt.ylabel(r'Transmission')
plt.savefig('./files_benzene/transmission.eps')
plt.close()

# Next, we wish to construct a model for the pi system. 
# This model should include all the other basis functions in terms of a self energy.
# The result is an effective "Hamiltonian" which is energy-dependent. 
plt.figure(1,(3,3))
H_eff , S_other = benz.get_effective_Hamiltonian(pzs, energies) 

np.set_printoptions(suppress=True, precision =2)
print "Effective Hamiltonian at Fermi energy"
print H_eff[:,:,len(energies)/2].real

# Plot energy dependence of effective Hamiltonian 
plt.plot(energies, H_eff[0,0,:], label = r'$\tilde{H}_{11}(E)$', c = 'b')
plt.plot(energies, H_eff[0,1,:], label = r'$\tilde{H}_{12}(E)$', c = 'r')
plt.plot(energies, H_eff[0,2,:], label = r'$\tilde{H}_{13}(E)$', c = 'g')
plt.plot(energies, H_eff[0,3,:], label = r'$\tilde{H}_{14}(E)$', c = 'k')

plt.legend(loc ='lower right')
plt.xlabel(r'$E-E_F$ (eV)')
plt.ylabel(r'Matrix element')

plt.savefig('./files_benzene/H_eff.png')
plt.close()
