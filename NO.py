import numpy as np

def lowdin(s_mm):
	"""Return S^-1/2"""
	eig, rot_mm = np.linalg.eig(s_mm) #find eigenvalues overlap matrix
	eig = np.abs(eig)
	rot_mm = np.dot(rot_mm / np.sqrt(eig), rot_mm.T.conj()) #rotation matrix S^-1/2. 
	return rot_mm

def FermiDirac(E,T):
	kb=8.6173324*10**(-5) #eV/K
	return 1/(1+np.exp((E)/(kb*T)))

class NaturalOrbitals:
	"""Natural Orbitals class.

	A post-processing tool to find different types of natural orbitals.
	The implementation is after Reed, A.; Weinstock, R.; Weinhold, F. J. Chem. Phys. 1985, 83 (1985), 735. """

	def __init__(self, h=None, s=None,D=None, ne=None, mol=None, basis_dic=None, model =None, core_and_valence=False):
		"""Required arguments: 
			'h':		Hamiltonian/Kohn-Sham operator in atom-centered basis.
			's':		Overlap matrix in atom-centered basis. 
			'ne':		Number of electron pairs. It is assumed that the first 'ne' molecular eigenfunctions of h have occupations 2.0. 
			'mol':		Atoms object from e.g. Atomistic Simulation Environment (ASE).  
			'basis_dic':Dictionary which links atom index to list of basis function indices. Basis function 
			indices associated with each atom must be consequtive.  
			'model':	String. If different from None, an Atoms object will be saved. The ordering of the atoms in this object will have the same ordering as the natural hybrids
			 obtained from get_natural_bonding_orbitals(). 
			'core_and_valence': If only valence electrons in calculation: set this to True.  """
		print "####################\nInitializing NaturalOrbitals class\n####################"
		self.h = h
		self.s = s
		self.basis_dic= basis_dic
		self.dim = h.shape[0]
		self.ne = ne
		self.mol = mol
		self.NAOs = False
		self.NHOs = False
		self.model_name = model
		self.core_and_valence =core_and_valence
		if model!=None:
			from ase import Atoms
			self.model = Atoms() 
		self.transmission_uptodate= False 
		self.currents_uptodate = False


	def get_NAO(self, threshold = 1.90):
		""" 'threshold': Threshold value for occupation of lone-pair/bound state. """
		from copy import deepcopy
		if self.model_name!='None':
			from ase import Atom
		print "\nGetting NAOs... \n"
		#normalize basis
		if np.abs((np.diagonal(self.s)-np.ones(self.dim))).max()>10**(-6):
			print "Overlap matrix is not orthogonal. Normalizing basis. "
			r = np.identity(self.dim)
			norms_n = self.s.diagonal()
			r /= np.sqrt(norms_n)
			self.h = np.dot(r.T.conj(), np.dot(self.h, r))
			self.s = np.dot(r.T.conj(), np.dot(self.s, r))

		##### construct D and P ###### 
		#D is C \lambda Cdagger, where C is the matrix of coefficients of the
		# natural orbitals of the wave functions (expanded in the basis {chi_k} 
		#and \lambda is the diagonal matrix of the occupation numbers of the natural orbitals)
		#find natural orbitals of full H by transforming to lowdin basis
		rot_lowdin  = lowdin(self.s) 
		Hlow = np.dot(rot_lowdin.T.conj(), np.dot(self.h, rot_lowdin))
		Slow = np.dot(rot_lowdin.T.conj(), np.dot(self.s, rot_lowdin))
		eig_v, eig = np.linalg.eig(Hlow)
		sort_list = eig_v.real.argsort()
		eig_v = eig_v[sort_list]
		eig = eig.take(sort_list, axis=1)
		norms_n = np.dot(np.dot(eig.T.conj(), Slow), eig).diagonal()
		eig /= np.sqrt(norms_n)
		C = np.dot(rot_lowdin, eig)
		print "Number of electron pairs:", self.ne
		self.HOMO = eig_v[self.ne-1]
		self.LUMO = eig_v[self.ne]
		self.E_F = 0.5*(eig_v[self.ne-1]+eig_v[self.ne]) # easy way to find midway point between HOMO and LUMO
		Lambda = np.zeros((self.dim, self.dim), dtype = 'complex')
		for n in range(self.ne):
			Lambda[n,n] = 2.0
		D = np.dot(C, np.dot(Lambda, C.T.conj()))
		self.D = D
		# P = SDS
		self.P = np.dot(self.s, np.dot(D, self.s))

		####### formation of PNAOs #######
		#make dictionary of saved indices on each atom
		d = {}
		for n in range(len(self.mol)):
			d[n] = 0
		PNAOs = np.identity(self.dim, dtype ='complex')
		self.PNHO = np.identity(self.dim, dtype ='complex')
		occ_PNAOs = np.identity(self.dim)
		self.PNHO_indices = []
		self.NMB_indices = []
		self.NRB_indices = []
		print "Getting PNAOs..."
		for n, atom in enumerate(self.mol):
		#	print "atom ", n, atom.symbol
			bfs =self.basis_dic[n]
			if atom.symbol in ['H']:
				self.NMB_indices.append(bfs[0])
				self.NRB_indices+=bfs[1:]
			if atom.symbol in ['N', 'C','O']:	
				if self.core_and_valence==True:
					self.NMB_indices+= bfs[:5]
					self.NRB_indices+=bfs[5:]
				else:
					self.NMB_indices+= bfs[:4]
					self.NRB_indices+=bfs[4:]
			if atom.symbol in ['S', 'Si', 'P']:	
				if self.core_and_valence==True:
					self.NMB_indices+= bfs[:9]
					self.NRB_indices+=bfs[9:]
				else:
					self.NMB_indices+= bfs[:4]
					self.NRB_indices+=bfs[4:]
			if atom.symbol in ['Au']:	
				self.NMB_indices+= bfs[:6]
				self.NRB_indices+=bfs[6:]

			p_sub = self.P.take(bfs,axis=0).take(bfs, axis=1)
			s_sub = self.s.take(bfs,axis=0).take(bfs, axis=1)	
			## lowdin orthogonalize
			rlow = lowdin(s_sub) #S^(-1/2)
			p_low = np.dot(rlow.T.conj(), np.dot(p_sub, rlow))
			s_low = np.dot(rlow.T.conj(), np.dot(s_sub, rlow))
			## solve generalized eigenvalue problem
			w_A_L, n_A_L = np.linalg.eig(p_low)
			sort_list = w_A_L.real.argsort()[::-1]
			w_A_L = w_A_L[sort_list]
			n_A_L = n_A_L.take(sort_list, axis=1)
			norms_n = np.dot(np.dot(n_A_L.T.conj(), s_low), n_A_L).diagonal()
			n_A_L /= np.sqrt(norms_n) #normalize lowdin basis
			n_A = np.dot(rlow, n_A_L) ### transform to local basis
			PNAOs[bfs[0]:bfs[-1]+1, bfs[0]:bfs[-1]+1] = n_A
			for n, bf in enumerate(bfs):
				occ_PNAOs[bf, bf] = w_A_L[n].real
		norms_n = np.dot(np.dot(PNAOs.T.conj(), self.s), PNAOs).diagonal()
		PNAOs /= np.sqrt(norms_n) # ensure normalization
		self.PNAOs = PNAOs


		######### Occupancy-weighted interatomic orthogonalization within NMB #####
		s_PNAO =  np.dot(self.PNAOs.T.conj(), np.dot(self.s, self.PNAOs))
		NAOs = deepcopy(self.PNAOs) # copy and work from PNAOs
		bfs = self.NMB_indices
		W_sub = occ_PNAOs.take(bfs,axis=0).take(bfs, axis=1)
		s_sub = s_PNAO.take(bfs,axis=0).take(bfs, axis=1)
		O_w = np.dot(W_sub, lowdin(np.dot(W_sub, np.dot(s_sub,W_sub.T.conj()))))
		test = NAOs.take(range(self.dim),axis=0).take(bfs, axis=1)
		for n, bf in enumerate(bfs):
			NAOs[:,bf] = np.dot(O_w[:,n],test.T.conj())	


		#### Gram-Schmidt orthogonalize unoccupied (NRB) onto occupied (NMB) ###
		print "\nIndices of Natural minimal basis (NMB):", self.NMB_indices
		print "Indices of Natural Rydberg basis (NRB)\n", self.NRB_indices
		s_PNAO = np.dot(NAOs.T.conj(), np.dot(self.s, NAOs))
		for n in self.NRB_indices:
			v = NAOs[:,n] # unocuppied
			for m in self.NMB_indices:#+ortho:# subtract projection onto orthogonalized vectors
				#		quit()
				v-=s_PNAO[n,m]*NAOs[:,m]
			norm = np.dot(np.dot(v.T.conj(), self.s), v)
			v/= np.sqrt(norm)
			NAOs[:,n] = v

		norms_n = np.dot(np.dot(NAOs.T.conj(), self.s), NAOs).diagonal()
		NAOs /= np.sqrt(norms_n) # ensure normalization

		####### Restoration of natural character of the NRB ####
		P_restore = np.dot(NAOs.T.conj(), np.dot(self.P, NAOs))
		s_restore = np.dot(NAOs.T.conj(), np.dot(self.s, NAOs))

		for n, atom in enumerate(self.mol):
			if len(self.NRB_indices)==0:
				break
#			print "atom ", n, atom.symbol
			bfs = self.basis_dic[n]
			bfs = list(set(bfs).intersection(self.NRB_indices)) 
			p_sub = P_restore.take(bfs,axis=0).take(bfs, axis=1)
			s_sub = s_restore.take(bfs,axis=0).take(bfs, axis=1)	
			## lowdin orthogonalize
			rlow = lowdin(s_sub) #S^(-1/2)
			#print rlow.shape
			p_low = np.dot(rlow.T.conj(), np.dot(p_sub, rlow))
			s_low = np.dot(rlow.T.conj(), np.dot(s_sub, rlow))
			## solve generalized eigenvalue problem
			w_A_L, n_A_L = np.linalg.eigh(p_low) #np.linalg.eigh ensures linearly independent eigenvectors. This is related to numerical problem indicated in section 4.a of Appendix. 
			sort_list = np.abs(w_A_L).argsort()[::-1]
			w_A_L = w_A_L[sort_list].real
			n_A_L = n_A_L.take(sort_list, axis=1)
			norms_n = np.dot(np.dot(n_A_L.T.conj(), s_low), n_A_L).diagonal()
			n_A_L /= np.sqrt(norms_n) #normalize lowdin basis
			n_A = np.dot(rlow, n_A_L) ### transform to basis prior to this loop
			test = NAOs.take(range(self.dim),axis=0).take(bfs, axis=1)
			for m, bf in enumerate(bfs):
			 	occ_PNAOs[bf, bf] = np.abs(w_A_L[m].real)
				if occ_PNAOs[bf, bf]<10**(-8): #
					occ_PNAOs[bf, bf] = 10**(-8)#resolve numerical problem as indicated in section 4.a of Appendix
				NAOs[:,bf] = np.dot(n_A[:,m],test.T.conj())	

		norms_n = np.dot(np.dot(NAOs.T.conj(), self.s), NAOs).diagonal()
		NAOs /= np.sqrt(norms_n) # ensure normalization

		####### occupancy weighted orthogonalization of NRB ######
		if len(self.NRB_indices)>0:
			self.s_PNAO =  np.dot(NAOs.T.conj(), np.dot(self.s, NAOs))
			bfs =self.NRB_indices
			W_sub = occ_PNAOs.take(bfs,axis=0).take(bfs, axis=1)
			s_sub = self.s_PNAO.take(bfs,axis=0).take(bfs, axis=1)
			O_w = np.dot(W_sub, lowdin(np.dot(W_sub, np.dot(s_sub,W_sub.T.conj()))))
			norms_n = np.dot(np.dot(O_w.T.conj(), s_sub), O_w).diagonal()
			O_w /= np.sqrt(norms_n) # ensure normalization
			test = NAOs.take(range(self.dim),axis=0).take(bfs, axis=1)
#			ok = []
#			problem = []
			for n, bf in enumerate(bfs): #express in local basis
#				print O_w[:,n].max(), bf
				NAOs[:,bf] = np.dot(O_w[:,n],test.T.conj())	
		norms_n = np.dot(np.dot(NAOs.T.conj(), self.s), NAOs).diagonal()
		NAOs /= np.sqrt(norms_n) # ensure normalization

		####### Formation of final NAOs ########

		P_test = np.dot(NAOs.T.conj(), np.dot(self.P, NAOs))
		s_test = np.dot(NAOs.T.conj(), np.dot(self.s, NAOs))

		for n, atom in enumerate(self.mol):
			bfs =self.basis_dic[n]
			p_sub = P_test.take(bfs,axis=0).take(bfs, axis=1)
			s_sub = s_test.take(bfs,axis=0).take(bfs, axis=1)	
			## solve eigenvalue problem
			w_A, n_A = np.linalg.eigh(p_sub)
			sort_list = w_A.real.argsort()[::-1]
			w_A = w_A[sort_list]
			n_A = n_A.take(sort_list, axis=1)
			norms_n = np.dot(n_A.T.conj(), np.dot(s_sub, n_A)).diagonal()
			n_A /= np.sqrt(norms_n) #normalize basis
			NAOs[:,bfs[0]:bfs[-1]+1] = np.dot(NAOs[:,bfs[0]:bfs[-1]+1], n_A)
			norms_n = np.dot(np.dot(NAOs.T.conj(), self.s), NAOs).diagonal()
			NAOs /= np.sqrt(norms_n) # ensure normalization
		#	# check for bound states
			#print atom.symbol, w_A.real
			for k, n_i in enumerate(w_A.real):
				if n_i<threshold:
					break
				if self.model_name!=None:
					self.model+= Atom('X', position=atom.position)
		#		print "found bound state", n_i
				self.PNHO[:, bfs[k]] = NAOs[:, bfs[k]]
				d[n] = d[n]+1
				self.PNHO_indices.append(bfs[k])
		self.d = d
		print "Bound states/lone-pairs (atom, number of bound states): ", self.d
		self.NAOs = NAOs[:,:]
		self.s_NAO = np.dot(self.NAOs.T.conj(), np.dot(self.s, self.NAOs))
		self.h_NAO = np.dot(self.NAOs.T.conj(), np.dot(self.h, self.NAOs))
		self.P_NAO = np.dot(self.NAOs.T.conj(), np.dot(self.P, self.NAOs))

		c = NAOs.take(self.NMB_indices+self.NRB_indices, axis=1)
		test = np.dot(c.T.conj(),np.dot(self.h, c))
		print "\nElectrons in occupied NAOs", np.trace(self.P_NAO.take(self.NMB_indices, axis=0).take(self.NMB_indices,axis=1)).real
		print "Total electrons input to calculator:", self.ne*2
		print "Total electrons calculated: :", np.diagonal(self.P_NAO).sum().real
		print "Percentage of electrons represented by NMB: %.8f \n" %(np.trace(self.P_NAO.take(self.NMB_indices, axis=0).take(self.NMB_indices,axis=1)).real/(2*self.ne)*100)
		return self.NAOs

	def get_natural_charges(self):
		if type(self.NAOs)==type(False):
			self.get_NAO()
		self.charges = np.zeros(len(self.mol))
		for n, atom in enumerate(self.mol):
			bfs = self.basis_dic[n]
			self.charges[n] = np.diagonal(self.P_NAO.real).take(bfs,axis=0).sum()
		return self.charges

	def get_NAO_occupation(self, NAO_indices):
		if type(self.NAOs)==type(False):
			self.get_NAO()
		return np.trace(self.P_NAO.take(NAO_indices, axis=0).take(NAO_indices,axis=1)).real

	def get_single_occupancy_NAO_indices(self, threshold=0.01):
		if type(self.NAOs)==type(False):
			self.get_NAO()
		occ = np.abs(np.diagonal(self.P_NAO)-np.ones(self.dim))
		indices = []
		for n, o in enumerate(occ):
			if o<threshold:
				indices.append(n)
		return indices

	def get_natural_bonding_orbitals(self, threshold=1.90):
		"""Implementation after Foster, J. P.; Weinhold, F. J. Am. Chem. Soc. 1980, 102 (22), 7211. """
		if self.model_name!=None:
			from ase import Atom
		if type(self.NAOs)==type(False):
			self.get_NAO()
		if type(self.NHOs)!=type(False):
			return self.NHOs
		###### deplete lone-pairs and bound states #####

		for n, atom in enumerate(self.mol): 
			bfs = self.basis_dic[n]
			p_sub = self.P_NAO.take(bfs,axis=0).take(bfs, axis=1)
			#deplete lone pairs and bound states
		#		print "Depleting", d[n]
			for l in range(self.d[n]):
#				print "l", l
				p_sub-=self.P_NAO[bfs[l],bfs[l] ]*np.dot(self.PNHO.take(bfs,axis=0).take([l], axis=1), self.PNHO.take(bfs,axis=0).take([l], axis=1).T.conj())
			self.P_NAO[bfs[0]:bfs[-1]+1,bfs[0]:bfs[-1]+1]= p_sub

		###### make hybrids #######
		print "Looking for electron pairs in bonds..."
		for n, atom in enumerate(self.mol): 
			#loop over pairs"
			for m, neigh in enumerate(self.mol):
				if n==m:
		#			print n, m, "break"
					continue
		#		print n, atom.symbol, m, neigh.symbol
				# find eigenvectors with high occuptions
				bfs1 = self.basis_dic[n]
				bfs2 = self.basis_dic[m]
				bfs = bfs1+bfs2
				p_sub = self.P_NAO.take(bfs,axis=0).take(bfs, axis=1)
				n_A, c_A = np.linalg.eig(p_sub)
				sort_list = n_A.real.argsort()[::-1]
				n_A = n_A[sort_list]
				c_A = c_A.take(sort_list, axis=1)
				for k, n_i in enumerate(n_A):
		#			print "Trying index", k,"bfs:",  bfs[k], "occ.:", n_i
					if n_i<threshold:
						break
		#			print "found hybrid", n, atom.symbol, m, neigh.symbol
					if self.model!=False:
						self.model+= Atom('X',atom.position+0.25*(neigh.position-atom.position))
					self.PNHO[bfs1[0]:bfs1[-1]+1, bfs[self.d[n]]] =np.dot(self.NAOs[bfs1[0]:bfs1[-1]+1,bfs1[0]:bfs1[-1]+1], c_A[:len(bfs1), k]) #rotate to local basis
					self.PNHO_indices.append(bfs1[self.d[n]])
					self.d[n] = self.d[n]+1
		if self.model!=False:
			from ase.io import write
			write('%s.traj'%self.model_name, self.model)

		norms_n = np.dot(np.dot(self.PNHO.T.conj(), self.s), self.PNHO).diagonal()
		self.PNHO /= np.sqrt(norms_n)
		print "Bound states and electron pairs in bonds: (atom, number of bound states) ", self.d

		#plot_basis( PNHO,mol, PNHO_indices , folder_name='PNHO', vacuum=3.0)
		#print PNHO_indices
		#quit()

		##### Lowdin orthogonalize #######
		print "Lowdin orthogonalizing NHO..."
		rot_low = lowdin(np.dot(np.dot(self.PNHO.take(self.PNHO_indices, axis=1).T.conj(), self.s), self.PNHO.take(self.PNHO_indices, axis=1)))
		self.NHOs = np.dot( self.PNHO.take(self.PNHO_indices, axis=1), rot_low)

		print "Finding unocuppied basis functions..."
		##### Gram-Schmidt orthogonalize #####
		rot = np.identity(self.dim, dtype ='complex')
		rot[:,:len(self.PNHO_indices)] = self.NHOs[:,:]
		for n in range(len(self.PNHO_indices),self.dim):
#			print "n", n
			v = rot[:,n]
			norm = np.dot(np.dot(v.T.conj(), self.s), v)
			v/= np.sqrt(norm)
			for m in range(n):# subtract projection onto orthogonalized vectors
				proj = np.dot(v.T.conj(), np.dot(self.s, rot[:,m]))
				v=v-proj*rot[:,m]
				norm = np.dot(np.dot(v.T.conj(), self.s), v)
				v/= np.sqrt(norm)
			rot[:,n] = v
		norms_n = np.dot(np.dot(rot.T.conj(), self.s), rot).diagonal()
		rot /= np.sqrt(norms_n)
		self.NHOs = rot
		self.h_NHO = np.dot(rot.T.conj(), np.dot(self.h, rot))
		self.s_NHO = np.dot(rot.T.conj(), np.dot(self.s, rot))
		print "Found natural hybrid orbitals. "
		return self.NHOs

	def get_NHO_hamiltonian_and_overlap(self):
		if type(self.NAOs)==type(False):
			self.get_NAO()
		if type(self.NHOs)==type(False):
			self.get_natural_bonding_orbitals()
		H = np.dot(self.NHOs.T.conj(), np.dot(self.h,self.NHOs))
		S = np.dot(self.NHOs.T.conj(), np.dot(self.s,self.NHOs)).real
		return H,S

	def get_NAO_retarded_Greens_function(self, energies):
		if type(self.NAOs)==type(False):
			self.get_NAO()
		G_NAO = np.zeros( (  self.dim, self.dim, len(energies) ), dtype = 'complex')
		for e, energy in enumerate(energies):
#			print e, energy
			G_NAO[:, :,e] = np.linalg.inv(self.s_NAO*energy - self.h_NAO)
		self.G_NAO = G_NAO
		return self.G_NAO

	def get_transmission(self, energies, SigmaL, SigmaR,h2=None, s2=None, eta=10**(-6)):
		if h2==None:
			if type(self.NAOs)==type(False):
				self.get_NAO()
			h2 = self.h_NAO
			s2 = self.s_NAO	
		self.energies=energies
		t = np.zeros(len(energies))
		GammaL = -2*SigmaL.imag.astype(np.dtype('complex'))
		GammaR = -2*SigmaR.imag.astype(np.dtype('complex'))
		print "Calculating transmission..."
		for e, energy in enumerate(energies):
#			print e, energy
			Gr = np.linalg.inv(-SigmaL-SigmaR + s2*(energy+eta) - h2)
			t[e] = np.dot(GammaL, np.dot(Gr, np.dot(GammaR, Gr.T.conj()))).real.trace()
		print "Calculated transmission."
		self.t = t
		self.transmission_uptodate= True
		return t
	def get_currents(self,Vs, T):
		"Return currents in units of 2e^2/h. "
		from scipy.integrate import simps
		if self.transmission_uptodate==False:
			print "You must calculate transmission first!"
			quit()
		Is = np.zeros(len(Vs))
		prefac = np.zeros(self.energies.shape)
		for v, V in enumerate(Vs):
		 	for e, energy in enumerate(self.energies):
		 		prefac[e] = FermiDirac(energy-self.E_F-V/2.0, T)-FermiDirac(energy-self.E_F+V/2.0, T)
			integrand = np.multiply(prefac, self.t) 	
			Is[v] = simps(integrand,self.energies)        
		self.Is = Is
		self.Vs = Vs
		self.currents_uptodate = True
		return Is

	def get_conductances(self):
		if self.currents_uptodate==False:
			print "You must calculate currents first!"
		Gs = np.zeros(self.Vs.shape)
		for v, V in enumerate(self.Vs):
			Gs[v] = self.Is[v]/V
		self.Gs = Gs
		return Gs


	def get_effective_Hamiltonian(self, indices, energies,h2=None, s2=None):
		""" 
		Function to calculate energy-dependent effective Hamiltonian. 
		The Green's function of the effective Hamiltonian reproduces the
		corresponding Green's function elements of the full Hamiltonian. 
		"""
		if h2==None:
			if type(self.NAOs)==type(False):
				self.get_NAO()
			h2 = self.h_NAO
			s2 = self.s_NAO			
		print "Getting effective Hamiltonian"
		h = h2.take(indices, axis=0).take(indices, axis=1)
		others = list(set(range(h2.shape[0])) - set(indices))
		s_other = s2.take(others,axis=0).take(others, axis=1) #must be complex
		h_other = h2.take(others, axis=0).take(others, axis=1)
		g_other = np.zeros((h_other.shape[0], h_other.shape[0], len(energies)), dtype = 'complex')
		h_tc = h2.take(others, axis=0).take(indices,axis=1)
		for n, e in enumerate(energies):
			g_other[:,:,n]= np.linalg.inv(e*s_other-h_other)
		S_other = np.zeros((len(indices),len(indices),len(energies)), dtype = 'complex')
		H_eff = np.zeros((len(indices),len(indices),len(energies)), dtype = 'complex')
		for n, e in enumerate(energies):
			S_other[:,:,n] = np.dot(h_tc.T.conj(), np.dot(g_other[:,:,n], h_tc))
			H_eff[:,:,n] =  h + S_other[:,:,n]
		return H_eff, S_other



def plot_basis(r, atoms, ns, basis = 'dzp', folder_name='./basis', vacuum=3.0, h = 0.20):
    """
    r: coefficients of atomcentered basis functions
    atoms: Atoms-object 
    ns: indices of bfs functions to plot. 

    """
    from gpaw import GPAW
    from gpaw import setup_paths
    from os.path import exists
    from subprocess import call
    from numpy import ascontiguousarray as asc
    from ase.io import write

    if exists(folder_name)==False:
        print "making folder for basis functions"
        call('mkdir %s'%folder_name, shell=True)

    symbols = atoms.get_chemical_symbols()
    ns=np.array(ns)
    atoms.center(vacuum=vacuum)
    calc = GPAW(h = h, mode='lcao',
        basis= basis,
        txt='basis.txt')
    atoms.set_calculator(calc)
    calc.initialize(atoms)
    calc.set_positions(atoms)
    c_fo_xi = asc(r.real.T)#coefficients
    phi_xG = calc.wfs.basis_functions.gd.zeros(len(c_fo_xi))
    calc.wfs.basis_functions.lcao_to_grid(c_fo_xi, phi_xG, -1)
    for n, phi in zip(ns, phi_xG.take(ns, axis=0)):
        print "writing %d of %d" %(n, len(ns)), 
        write('%s/%d.cube' % (folder_name,n), atoms, data=phi)
    summ = np.zeros(phi_xG[0,:,:,:].shape)
    print "sum", summ.shape
    for n, phi in zip(ns, phi_xG.take(ns, axis=0)):
        summ+= phi
    write('%s/sum.cube'%folder_name, atoms, data=summ)


def rwb(mini, maxi, value):
    avg = mini+(maxi-mini)/2.0
    #print "avg", avg
    #print value
    if value>=avg:
        red = 255
        green = int(-255.0/(maxi-avg)*(value-maxi))
        blue = green
    if value<avg:
        blue =255
        green =int(255/(avg-mini)*(value-mini))
        red = green
#   print red, green, blue
    return red, green, blue

def get_basis_STOnG(symbols, indices, basis):
    nZeta = len(basis.split('-')[1].split('+')[0].split('G')[0])
    #if 'd' in basis.split('-')[1].split('+')[0].split('G')
    #print basis.split('-')[1].split('+')[0].split('G')[0]
    pol_H =3*basis.count(',')
    pold =5*basis.count(',')
    diff_sp =0
    diff_s = 0
    if basis.count('+')==1:
        diff_sp = 4
    print basis
    if basis.count('+')==2:
        diff_sp = 4
        diff_s = 1
    bfs = []
    count =0
    # print "nZeta", nZeta
    # print "pol_H", pol_H
    # print "d-polarization", pold
    # print "diffuse sp", diff_sp
    # print "diffuse s", diff_s
    for n, atom in enumerate(symbols):
        if atom in ['H']:
            if n in indices:
                bfs+=range(count, count+nZeta+diff_s+pol_H)
            count+=nZeta+diff_s+pol_H
        if atom in ['C', 'N', 'O', 'B', 'F']:
            if n in indices:
                bfs+= range(count, count+1+4*nZeta+pold+diff_sp)
            count+=1+4*nZeta+pold+diff_sp
        if atom in ['S', 'P']:
            if n in indices:
                bfs+= range(count, count+1+4+4*nZeta+pold+diff_sp)
            count+=1+4+4*nZeta+pold+diff_sp
    return bfs



def jmol_effective_Hamiltonian(mol, Heff, outfile = 'H_eff', cutoff=0.5):
	###### plot direct NHO Hamiltonian ######
	mol.center()
	from ase.io import write
	maxi= np.abs(Heff).max()#H_eff.min()
	mini = -maxi
	dummy = mol[:]
	Heff = Heff.real
	dummy.translate(np.array([[5.0, 0.0, 0.0]]))
	write('%s.xyz'%outfile, dummy)
	symbols = mol.get_chemical_symbols()
	scaling_factor = 8.0
	f = open('%s.dat'%outfile, 'w')
	f.write('#test \n')
	f.write('wireframe off \n')# background white \n')
	for n, orbital in enumerate(mol):
		onsite =Heff[n, n] #, len(energies)/2]
	#	print onsite
		diam = 0.5#np.abs(onsite)/3.0
		print "diam", diam
		red, green, blue = rwb(mini, maxi, onsite)
		f.write('draw sphere%d diameter %f {%f %f %f} color {%d %d %d }\n'%(orbital.index, diam, orbital.position[0],orbital.position[1], orbital.position[2], red , green , blue)) 
		for m, neigh in enumerate(mol):
			#if neigh.symbol =='Au' or neigh.symbol =='H' or neigh.index ==atom.index:	
			#		continue
			coupling=Heff[m,n]#, len(energies)/2 ]
			r,g,b = rwb(mini, maxi, np.abs(coupling))
	#		if np.abs(coupling)>2.3:
	#			print h[m,n]
	#		if n==6 and m==7:
	#			print "<6 H 7 >=", coupling
	#		if n==6 and m==0:
	#			print "<6 H 7 >=", coupling

		 	if np.abs(coupling)>cutoff:					
					f.write('draw arrow%d_%d arrow {%f %f %f} {%f %f %f} diameter %.6f color {%f %f %f}\n'%(orbital.index,neigh.index,  orbital.position[0],orbital.position[1], orbital.position[2], neigh.position[0],neigh.position[1], neigh.position[2], np.abs(coupling)/(4*6.0), r,g,b)) 



