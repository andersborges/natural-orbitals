# From computational experiment to quantifiable understanding

This implementation of natural orbitals can be used to e.g.:

1. Extract parameters for Huckel-type/tight-binding models. 
2. Calculate the electrical conductance of a molecule from a gas phase calculation. 

![alt text](pic.png)


# What is it?
natural-orbitals is a Python class written to calculate different types of natural orbitals based on the output of electronic structure codes. It can also calculate effective Hamiltonians and calculate the electrical conductance on the basis of the Landauer formula. 

# Usage
Run an example file with python: `python benzene_transport_example.py` or go through the tutorial belown. 

# How does it work?
Supply the class with input extracted from DFT or HF calculations. Input from [GPAW](https://wiki.fysik.dtu.dk/gpaw/) is used in the example but output from Gaussian09 is also possible. 

The implementation of natural atomic orbitals (NAO) is based on: 

[Reed, A.; Weinstock, R.; Weinhold, F. _J. Chem. Phys._ 1985, 83 (1985), 735.](http://dx.doi.org/10.1063/1.449486)

Note: Spherical symmetry is not implemented. As a result, the NAOs do not neccesarily resemble s, p or d-orbitals. 

The implementation of natural hybrid orbitals is based on: 

[Foster, J. P.; Weinhold, F. _J. Am. Chem. Soc._ 1980, 102 (22), 7211.](http://dx.doi.org/10.1021/ja00544a007)

# Tutorial
This tutorial will go through some of the features of the class. As an example we will look at the benzene molecule. 

## Initializing class
The class can be initialized with the following parameters. 

```python 
		"""Arguments:  
			'h':		(N,N) ndarray
						Hamiltonian/Kohn-Sham operator in atom-centered basis.
			's':		(N,N) ndarray
						Overlap matrix in atom-centered basis. 
			'ne':		integer int
						Number of electron pairs. It is assumed that the first 'ne' molecular eigenfunctions of h have occupations 2.0. 
			'mol':		
						Atoms object from e.g. Atomistic Simulation Environment (ASE).
			'basis_dic':	dict
					Dictionary which links atom index to list of basis function indices. Basis function indices associated with each atom must be consequtive.  
			'model':	
						{None, str}, optional
						If different from None, an Atoms object will be saved. The ordering of the atoms in this object will have the same ordering as the natural hybrids
			 obtained from get_natural_bonding_orbitals(). 
			'core_and_valence': {None, bool)
						If input is for an all-electron calculation (like Gaussian): set this to True."""

```

## Obtaining input
In principle, input from any electronic structure code can be used.
#### Input from GPAW
The tool dump_hamiltonian_parallel() can be used to save the Hamiltonian and overlap matrix to disk. If calc is a converged GPAW calculator and mol is an atoms object:

```python
from gpaw.lcao.tools import dump_hamiltonian_parallel as dhp
calc.dump_hamiltonian_parallel('hs', mol)
```
This will generate a pickle file. The library basis_dic can be generated with e.g. :
```python
from gpaw.lcao.tools import get_bfi2
basis_dic = {}
for atom in mol:
	basis_dic[atom.index]=get_bfi2(mol.get_chemical_symbols(), 'dzp', [atom.index])
```

#### Input from Gaussian09
Obtaining input from Gaussian is slightly more tricky. One way to go about it is to save the read-write file (.rwf) by including the following lines in the beginning of the Gaussian input file (.com):
```
%rwf=mol.rwf
%save
```
The .rwf can then copied to the working directory. The script parse_g09.sh can then be used to extract the required input. The script takes the output file of the calculation as input. This script relies on the tool rwfdump.  Use tool directly from the terminal: 
```python parse_g09.py mol.log```

This will create a folder containing the necessary input. 

## Obtaining natural orbitals
After the class has been initialized, different kinds of natural orbitals can be obtained as linear combination of the initial local basis. 

By calling .get_NAO() you can get:
* natural molecular orbitals: (.NO)
* pre-natural atomic orbitals (.PNAO)
* natural atomic orbitals (.NAOs)

By calling .get_natural_bonding_orbitals() you can get:
* pre-natural hybrid orbitals (.PNHO)
* natural hybrid orbitals (.NHOs)

You can use the function plot_basis() to plot the orbitals obtained from GPAW calculations:
```python
plot_basis(NO.NAOs, mol, pzs, basis = 'dzp', folder_name='./files_benzene/pzs') 
```

## Obtaining matrix elements
You can obtain the Hamiltonian and overlap matrix in your desired basis by rotating them yourself. Here is an example:

```python
h_NHO = np.dot(NO.NHOs.T.conj(), np.dot(h, NO.NHOs))
```

## Calculating transport properties
Is implemented. Will write more later. 

