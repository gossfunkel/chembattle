from ursina import *
from chempy import Substance
from random import uniform
from random import randint
import numpy as np
import atoms as ats

# ''' --NOTES: 
# Coulomb's law: f = k*q1*q2/r**
#k = 0.000000008988 # Nm**/C**
# Lennard-Jones potential is worked out via:
# Van der Waals' attractive force: F a r**6
# Pauli's repulsive force: F a 1/x**12
# '''
# ''' ---from ursina api ref for CubicBezier
# Curves from Ursina are used by Entity when animating, like this:
# e = Entity()
# e.animate_y(1, curve=curve.in_expo)
# e2 = Entity(x=1.5)
# e2.animate_y(1, curve=curve.CubicBezier(0,.7,1,.3))
# '''
# This is the game's physics engine.
# 
# I am in the process of rewriting a basic Lennard-Jones, Spring Bonding and Electric molecular dynamics sim. 
# 
# Don't expect this code to run reliably. The rest of the game isn't built up to keep up with these changes yet. This physics engine 
# and the networking class will be the main core of the game engine for this project, and if I can get a working prototype without 
# wasting the whole summer on it, I think it might be worth keeping up. 
# '''

# global values for simulation - can be modified in action
n 			= 0 	# number of particles in sim
D 			= 3 # number of spacial dimensions
setdt 		= 2.0E-5
LL 			= 25 	# scale of the simulation space in angstroms
L 			= np.zeros([D])+LL
Temp0 		= 200 	# maintained temperature in K - TODO: CHANGE FOR PRESSURE (walls)
#dragstate 	= False # bool value to pause the sim while dragging a molecule around

# constants
Kr  = 148000.0  #spring potential is usually defined as U = (k/2)(r-r_0)^2. can include /2
Kth = 35300.0  #same explanation as Kr but with bending energy
kb  = 0.8314459920816467 # Boltzmann constant in A^2 D ps^-2 K^-1
NA  = 6.0221409e+26 #Avogadros constant x 1000 (g->kg)
ech = 1.60217662E-19 #electron charge in coulombs
kc  = 8.9875517923E9*NA*1E30*ech*ech/1E24 #electrostatic constant in Daltons, electron charges, picosecond, angstrom units

# arrays of generated particles
#particles 		= []
atomPositions 	= [] 	# will be later initialised into a fast numpy array of vectors
atomVelocities 	= [] 	# will be later initialised into a fast numpy array of vectors
molecules 		= []	# list storing molecule objects
tp 				= [] 	# list tracking species of atom (element/isotope)
masses 			= [] 	# list tracking masses
#a 				= [] 	# array for use by update function to index arrays of molecules' components' accelerations
sigm 			= []	# array for population with sigma (balanced potential) LJP values for loaded atoms
epsi 			= []	# array for population with epsilon (potential well) LJP values for loaded atoms
#chrgs 			= [] 	# array for use by coul function to index all charges

# todo: properly outline boundaries
BC = False # True: periodic; False: reflective

def GetMolecules(player=0):
	return molecules
	# return a 1D array of all the molecules for a player
	# default to player 1 (local)
	#if (molecules.shape[player] > 0):
	#	return molecules[player, :]
	#else: return None

# smoothly transition a vector via linear interpolation
def SlideTo(pops, fpos, lerp):
	#target = Vec3([pops])
	#follow = Vec3([fpos])
	distance = np.linalg.norm(pops-fpos)
	if (distance > 0):
		direction = (pops - fpos) / distance
		min_step = max(0, distance - 100) #final int is max follow distance setting
		max_step = distance
		step_distance = min_step + (max_step - min_step) * lerp #LERP factor
		new_follow = Vec3(fpos + direction * (step_distance * time.dt))
	else: return(fpos)
	return(new_follow)

# CAUTION : OBJECT CAN BE INDEXED AS AN ARRAY-LIKE DATA TYPE THROUGH THE __array__() FUNCTION, SO BE CAREFUL HOW YOU CALL IT!
class Molecule(Entity):
	def __init__(self, nam, position, indx, chldatoms, tpis, velocity=np.zeros(D), bl=[0.95], th0=109.47*np.pi/180.0):
		if (len(chldatoms) == 0): return None # don't construct an empty molecule; what is that??? how would it even work?!
		self.n 			  = len(chldatoms)
		#self.size 		  = 4
		super().__init__(world_position=position)
		self.name 		  = nam
		self.visible_self = False # molecule isn't visible; atoms are visible
		#if self.n > 0:
		self.children 	  = chldatoms # the children of the molecule object are the atoms that constitute the molecule
		#print(self.children)
		self.indx 		  = indx # the index of the first child atom in the position and velocity arrays
		self.bl 		  = [bl] 	# 0.9611 for water -- bond lengths stored in the arrays stored in covbonds
		self.mm 		  = [] 	# array storing masses
		
		# construction lists
		covbonds		  = []	# array tracking covalent bonds between atoms; indexed by particle index
			# structure: [indx,bonded-to,bond-length,gamma(bond-strength)]
		chrg			  = [] 
		positions 		  = []
		velocities 		  = []
		angles		 	  = []	# dict tracking bond angles; indexed by particle index

		# iterate through and initialise values for all children (atoms in molecule)
		for i in range(self.n):
			# temporarily, this is set to extend the bond length array with duplicate values
			if i > len(bl)-1:
				bl.append(bl[i-1])
			if i > 0:
				# increment temporary position variable by bond length to give smoother load-in
				position = np.array(position) + np.array([1.0/bl[i]**3,1.0/bl[i]**3,1.0/bl[i]**3]) 
				# in future, this will be replaced with / followed by Maxwell-Boltzmann to match initial velocities to Temp
			else:
				position = np.array(position)
			positions.append(position)
			#print("adding atom at pos: ")
			#print(position)
			# initialise atom type, constants, absolute and relative positions, velocities, masses, charges, bonds, and angles
			child = self.children[i]
			#print(child)
			child.setIndx(self.indx+i)
			# the set of atoms in the arrays start at the molecule index; the molecule indexes the first atom in the array
			child.world_position = position
			child.velocity 		 = velocity
			#print("velocity:")
			#print(velocity[0])
			velocities.append(velocity)
			self.mm.append(child.mass)
			chrg.append(child.charge)
			# assign the atom's parentMol variable to the new molecule
			child.parentMol = self
			# initialise bonds
			if (self.n > 1):
				covbonds.append([bl[i],Kr])
				if (i > 0 and i < (self.n-1)):
					angles.append([self.children[i-1].indx,self.children[i].indx,self.children[i+1].indx,th0,Kth])

		# hacky temporary assigning of partial charge based on electronegativity
		# calculate partial charges using Allen electronegativity and formal charge contributions
		chrg = np.array(chrg)
		sumcharge = 0
		elnegal = []
		sumelecneg = 0
		for i in range(self.n):
			# count up the electronegativity of the atoms, and the formal charge on the molecule
			elnegal.append(self.children[i].electroNegAllen)
			sumelecneg += elnegal[i]
			sumcharge += chrg[i]
		# might use this? if i don't ditch all of this for a better method
		#sortedChrg = np.array(chrg)
		#sortedChrg.sort()
		#diffChrg = sortedChrg[-1] - sortedChrg[0]

		# find the difference between the largest and smallest electronegativities in the atom
		sortedElnegal = np.array(elnegal)
		sortedElnegal.sort()
		#print("sorted elnegal: ")
		#print(sortedElnegal)
		diffNeg = sortedElnegal[-1] - sortedElnegal[0]
		self.charge = sumcharge
		# go through children and calculate partial charge as a function of electronegativity
		for i in range(n):
			if (chrg[i] < 0): 	# assign sign based on charge of free atom
				diffNeg *= -1
			if sumcharge == 0:  # this molecule has a formal charge of 0
				chrg[i] = diffNeg * (elnegal[i] / sumelecneg) # latter term represents ratio of electronegativity in molecule
			else: # this molecule has an overall charge on it; it's an ion
				chrg[i] = sumcharge + diffNeg * (elnegal[i] / sumelecneg) # latter term represents ratio of electronegativity in molecule
			# AMBER uses e = sum(qiqj/rij) so could maybe rearrange. distance could define how much electrons can split charge
			# where at a distance beyond the maximum bonding distance the charge on an atom reaches the formal charge of a free atom
			# this might have to be saved to see if I can figure out a holistic approach to electrons
			self.children[i].charge = chrg[i] # set the charge value in the atom object to the partial charge
					
		#print(covbonds)
		# load lists into arrays
		print("adding the following to positions array:")
		print(np.array(positions))
		RescalePosArray(np.array(positions))
		RescaleVelArray([veloc for veloc in velocities])
		RescaleMassArray(self.mm)
		self.charges 		= chrg
		self.bl 			= np.array(bl)

		self.covbonds 	= np.array(covbonds)
		# bond angles form between two bonds, so are not required if there are less than that
		self.angles = angles # alternatively, np.array(angles.values())
		
		#print("molecule angles: ")
		#print(self.angles)

	def __array__():
		sumarray = np.array([self.children,self.mm,self.charges])
		return sumarray

	def update(self):
		#global atomPositions
		#global atomVelocities
		#global atomAccel
		for i in range(self.n):
			#print("index " + str(self.indx+i))
			print(atomPositions)
			self.children[i].world_position = SlideTo(atomPositions[self.indx+i, :],self.children[i].world_position,1)
			self.children[i].velocity 		= atomVelocities[self.indx+i, :]

	def positions(self):	return [self.children[i].world_position for i in range(self.n)]

	def velocities(self):	return [self.children[i].velocity for i in range(self.n)]

def reflectBC(r,v):
	newv = 1.0 * v
	newr = 1.0 * r
	print("r:")
	print(newr)
	for mol in molecules:
		for at in mol.children:
			for j in range(D):
				if newr[at.indx,j] < -L[j]:
					newr[at.indx,j]= -newr[at.indx,j]
					newv[at.indx,j]= abs(v[at.indx,j])
				if newr[at.indx,j] > L[j]:
					newr[at.indx,j]= 2.0*L[j]-newr[at.indx,j]
					newv[at.indx,j]= -abs(v[at.indx,j])
	return newr,newv

def RescalePosArray(*posits):
	global atomPositions
	global n
	# 	rescale atomPositions array
	if len(atomPositions) > 0:
		atmPos = atomPositions.tolist()
	else: atmPos = []
	# for each position array passed
	for pos in posits:
		#print(pos)
		# iterate through the positions in the array passed and append them to the atomPositions array
		atmPos.extend(pos)
	atomPositions = np.array(atmPos)
	n = len(atomPositions)
	#print(atomPositions)

def RescaleVelArray(velocs):
	global atomVelocities
	# 	rescale velocities array
	if len(atomVelocities) > 0:
		atmVel = atomVelocities.tolist()
	else: atmVel = []
	# for each velocity array passed
	for vel in velocs:	# iterate through the velocities in the array passed and append them to the atomVelocities array
		atmVel.extend(vel)
	atomVelocities = np.array(atmVel)

def RescaleMassArray(*newmas):
	global masses
	# 	rescale velocities array
	if len(masses) > 0:
		masss = masses.tolist()
	else: masss = []
	for mas in newmas:
		masss.extend(mas)
	masses = np.array(masss)

def RescaleSigmArray(*newsig):
	global sigm
	# 	rescale velocities array
	if len(sigm) > 0:
		si = sigm.tolist()
	else: si = []
	for sii in newsig:
		si.extend(sii)
	sigm = np.array(si)

def RescaleEpsiArray(*neweps):
	global epsi
	# 	rescale velocities array
	if len(epsi) > 0:
		ep = epsi.tolist()
	else: ep = []
	for epp in neweps:
		ep.extend(epp)
	epsi = np.array(ep)

def CreateMolecule(location, *mols, player=0, bl=[1], velocity=[np.zeros(D)]):
	global molecules
	global n

	atfact = ats.AtomFactory() # factory method for generating atoms

	# update location based on player info
	if (player == 0 and location == None): location = Vec3(18+randint(-5,5),5+randint(-1,1),8+randint(-5,5))
	elif (player == 1 and location == None): location = Vec3(32+randint(-5,5),5+randint(-1,1),92+randint(-5,5))

	# run the creation procedure for every molecule passed in *mols tuple
	# LOOP PER MOL PASSED TO CONSTRUCT
	for mol in mols:
		name = mol
		newatts = [] # instantiate list for new atoms in molecule
		newmol = Substance.from_formula(mol) # get a new chempy Substance object from the name
		#print("newmol = " + newmol.)
		comp = newmol.composition
		for elementindex,quant in comp.items(): # iterate through the members of the composition dictionary
			# print("key: " + str(elementindex) + ", value: " + str(quant))
			for i in range(quant): # loop for quantity
				# generate atoms of given type
				newt = atfact.createAtom(elementindex,player)
				newatts.append(newt) # add new atom to list
				print("created " + newt.name)
		# temporary messy shuffle to put 3-member molecules the right way round
		if len(newatts) > 2:
			tempAt = newatts[1]
			newatts[1] = newatts[2]
			newatts[2] = tempAt

		# initialise arrays for sigma, epsilon, and type index values for the new atoms in the molecule
		newsig 	= []
		neweps 	= []
		tpis 	= []
		
		# IS THE SIMULATION EMPTY?
		if n > 0:
			# if the simulation contains molecules, generate new LJ parameters using the Lorentz-Berthelot combining rules
			for i in range(len(newatts)):
				# 0 generate type name
				nam = mol + "_" + str(newatts[i].name)
				# 1 initialise new index
				tpindex  = len(tp)
				# 2 search typearray for atom type
				for typ in range(len(tp)):
					if tp[typ] == nam:
						tpindex = typ
				# type not found loaded: load new type, leave atom.unique = True
				if tpindex == len(tp):
					# 3 add type to array when loaded
					tp.append(nam)
					tpis.append(tpindex)
					newatts.tpindex = tpindex
					# 4 load simulation constants to arrays
					newsig.append(newatts[i].sig)
					neweps.append(newatts[i].eps)
					atsig = [newatts[i].sig]
					ateps = [newatts[i].eps]
					# si[0] and ep[0] should countain sigma and epsilon ii - the atom's LJ constants
					atsig.append([(newsig[i]+si[0])/2   for si in sigm]) # Lorentz combining rule: sigma combines using mean average. Analytically correct for hard spheres only
					ateps.append([(neweps[i]*ep[0])**.5 for ep in epsi]) # Berthelot combining rule: dipole interactions averaged with root-product
					# load newly calculated arrays for atom i as a dimension of the new arrays
					newsig[i] = atsig
					neweps[i] = ateps
					# update values of other mols in sim. The first member of every list in the sigm array contains the constant for each atom
					for si in sigm: si.tolist().append((si[0]+newatts[i].sig)/2) 	# Lorentz: (sig_ii+sig_jj)/2  		- mean
					for ep in epsi: ep.tolist().append((ep[0]*newatts[i].eps)**.5) 	# Berthelot: sqrt(eps_ii*eps_jj) 	- 3D pythagoras
				else: # type found in tp array: save type index to atom object for easy lookup and log that atom is only one of its type
					newatts[i].tpindex = tpindex
					tpis.append(tpindex)
					# tell engine that atom is not unique in the simulation
					newatts[i].setUnique(False)
		else: # if this is the first molecule added to the sim, initialise arrays to molecule values and leave atom.unique = True
			for i in range(len(newatts)):
				# 0 generate type name
				nam = mol + "_" + str(newatts[i].name)
				# 1 initialise new index
				tpindex  = len(tp)
				# 2 search typearray for atom type
				for typ in range(len(tp)):
					if tp[typ] == nam:
						tpindex = typ
				# type not found loaded: load new type, leave atom.unique = True
				if tpindex == len(tp):
					# 1 save indices
					newatts[i].tpindex = tpindex
					tp.append(nam)
					tpis.append(tpindex)
					# 2 load values - ARE THESE LOADED CORRECTLY?
					newsig.append([newatts[i].sig])
					neweps.append([newatts[i].eps])
				else: # type already loaded, copy index value
					newatts[i].tpindex = tpindex
					tpis.append(tpindex)
					# tell engine that atom is not unique in the simulation
					newatts[i].setUnique(False)
		# > load any new values into AMP arrays
		RescaleSigmArray(newsig)
		RescaleEpsiArray(neweps)

		# initiate a molecule with parameters passed
		mol = Molecule(name, location, len(molecules), newatts, tpis, velocity)
		#print(at.world_position)

		# add remaining values to AMP arrays and lists:
		molecules.append(mol)
		#print(mol.nam)
		#RescalePosArray(mol.positions())
		#RescaleVelArray(mol.velocities())
		# update number of atoms in simulation
		#n += mol.n DONE IN RESCALE FUNCTIONS

	# 	type, bond length & angle, sigma and epsilon value, mass, and charge arrays initialised in the molecule class
	# TODO: add to correct axis of arrays for which player is adding the molecule, default to player 1 (local) 

	# return the molecule object for any further assignments etc
	return newatts

#def ReactMols():
#	pass

# Gradient of Lennard-Jones potential on one body (negative force)
def dLJp(r,mol,atid,sigl,epsl,bdln):
	#atm = mol.children[atid]
	pvi = mol.children[atid].indx # position and vector array index of atom

	#print(mol.children[atid].name)
	#print(mol.children[atid].world_position)
	#print(mol.children[atid].velocity)
	#print(mol.children[atid].sig)
	#print(mol.children[atid].tpindex)
	sg = np.array([sigl])
	# if this is the only atom of its type, remove its sigma value to create the temporary array
	if mol.children[atid].isUnique(): sg=np.delete(sigl,0)
	# otherwise index 0 is used for interactions w/ others of its type

	ep = np.array([epsl])
	# remove epsilon value for self if it is unique
	if mol.children[atid].isUnique(): ep=np.delete(ep,0) 

	# format arrays to construct one with a value for every atom
	sixConst = ep*(sg**6)
	twelveConst = 2.0*ep*(sg**12)
	
	# calculate distance vector then absolute distance in each dimension
	drv=r-r[pvi] # distance vector for every point
	#print ("DRV:")print(drv)
	for j in range(mol.n): # remove position vectors of all atoms in the molecule from drv
		drv= np.delete(drv,mol.indx,0) # n.b. tge index of the removed value doesn't change between iterations due to the resizing of the array
	#print ("DRV after molecule removal:")print(drv)
	dr=[np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) for a in drv] #absolute distances of that lad

	# truncate to only calculate only for nearby (less than 4 fm away) atoms. May be replaced by partitioning
	if ((dr[i] > 4) for i in dr):
		dLJP=0.0
	else: # calculate dLJP (8-14 from 6-12)
		#print("ep: " + str(ep))
		#print("dr: " + str(dr))
		# TODO
		# 	figure out how the hell to look up the correct sigma and epsilon values for the target interactions
		for tpi in mol.tpis:
			# go through each type index in the simulation and run the calculations for atoms of this type
			r8   = sixConst[tpi] * (1.0/np.array(dr))**8 # IF THERE'S A DIVIDE BY ZERO, you're calculating particle self-interactions
			r14  = twelveConst[tpi] * (1.0/np.array(dr))**14
		r8v  = np.transpose(np.transpose(drv)*r8)
		r14v = np.transpose(np.transpose(drv)*r14)
		r8vs = np.sum(r8v,axis=0)
		r14vs= np.sum(r14v,axis=0)
		dLJP=24.0*(r14vs-r8vs)

	return dLJP

#gradient of bond length potential (negative force)
def dBEpot(r,mol,bnds):
	bps=np.zeros([mol.n,3])
	if mol.n > 1:
		for i in range(mol.n): #loop over particles in molecule
			#print("i: " + str(i) + ", n: " + str(mol.n))
			dr0=bnds[i][0]
			e0 =bnds[i][1]
			if (i == mol.n-1):
				dr=r[mol.children[i].indx]-r[mol.children[i-1].indx] # difference in position
			else: 
				dr=r[mol.children[i+1].indx]-r[mol.children[i].indx] # difference in position
			#print("dr from dbepot:")
			#print(dr)
			dr2=dr*dr 			# squared
			adr2=sum(dr2) 		# summed
			adr=np.sqrt(adr2) 	# rooted
			dBE=2.0*e0*(adr-dr0)*dr/adr 
			bps[i]+=dBE
	return bps

#gradient of bond angle potential (negative force)
def dBA(r,mol,angs):
	aps=np.zeros([mol.n,3])
	if mol.n > 1:
		for i in range(mol.n):
			for i in range(len(angs)): #loop over angles
				#i -= 1
				a1=int(angs[i][0])
				a2=int(angs[i][1])
				a3=int(angs[i][2])
				#print("dBA for " + mol.name)
				#print("a1: " + str(a1) + ", a2: " + str(a2) + ", a3: " + str(a3))
				childIndex = mol.children[i].indx
				if childIndex==a1 or childIndex==a2 or childIndex==a3:
					th00=angs[i][3] #equilibrium angle
					e0 =angs[i][4] #bending modulus
					if childIndex==a1 or childIndex==a2:
						r1=r[a1]-r[a2] #bond vector 1 (from middle atom to atom 1)
						r2=r[a3]-r[a2] #bond vector 2 (middle atom to atom 2)
					else:
						r1=r[a3]-r[a2] #bond vector 1 (from middle atom to atom 1)
						r2=r[a1]-r[a2] #bond vector 2 (middle atom to atom 2)
					print("r values in angle calculation: ra1:" + str(r[a1]) + ", ra2:" + str(r[a2]) + ", ra3:" + str(r[a3]))
					print("r values in angle calculation: r1:" + str(r1) + ", r2:" + str(r2))
					ar1=np.linalg.norm(r1) #lengths of bonds
					ar2=np.linalg.norm(r2)
					print("lengths of bonds in dBA: ar1: " + str(ar1) + ", ar2: " + str(ar2))
					dot=np.dot(r1,r2) #r1 dot r2
					print("dot product: " + str(dot))
					ndot=dot/np.linalg.norm(dot) #normalize dot product by vector lengths i.e. get the cos of angle
					print("normalised dot product in angle calculation: " + str(ndot))
					th=np.arccos(ndot) #bond angle, theta
					dUdth=-2.0*e0*(th-th00) #-dU/dtheta
					if a1==childIndex or a3==childIndex:
						numerator 	= (r2/(ar1*ar2))-(dot/(ar1**3 * ar2 * 2.0))
						denominator = np.sqrt(1.0 - ndot**2)
						dUdr 		= dUdth * numerator/denominator
						aps[i+1] 	+= dUdr
					if childIndex==a2:
						denominator	= np.sqrt(1.0 - ndot**2)
						#print (" denominator in angle calculation: " + str(denominator))
						n1 			= - (r2 + r1)
						n2			= dot * r1 / ar1**2
						n3			= dot * r2 / ar2**2
						numerator	= (n1+n2+n3)/(ar1*ar2)
						dUdr		= dUdth * numerator/denominator
						aps[i+1]	+= dUdr
						aps[i]		+= dUdr
						aps[i-1]	+= dUdr
	return aps

#derivative of coulomb potential (negative force)
def coul(r,molindex,i,q0,mol,mols):
	qs = []
	for ml in mols:
		if ml != mol:
			qs.extend(ml.charges)
		else: 
			# no charge interactions with own molecule
			#for c in ml.charges:
				#qs.append(np.zeros(3))
			modmolcharge = np.array(ml.charges)
			np.delete(modmolcharge,molindex)
			qs.extend(modmolcharge.tolist())
	#qs=1.0*np.array([ml.charges for ml in mols])
	qs = 1.0*np.array(qs)
	qs=np.delete(qs,i)
	#r = np.delete(r,i)
	drv=r-r[i] #distance in each dimension
	drv=np.delete(drv,i,0) #remove ith element (no self LJ interactions)
	dr=[np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) for a in drv] #absolute distance
	r3=q0*qs*kc*((1.0/np.array(dr))**3.0) # Coulomb's law
	FF =np.transpose(np.transpose(drv)*r3) # transpose the distance array, multiply by force, transpose back
	Fs=np.sum(FF,axis=0) # sum the value of axis 0 of the transposed & multiplied array
	#print(Fs)
	return Fs
	#return np.zeros(D) disable charge force

def rescaleT(v,T):
	KE = 0.5 * sum(sum(masses*np.transpose(v*v)))
	avKE = KE/n
	Tnow = (2.0/3)*avKE/kb
	lam=np.sqrt(T/Tnow)
	lam=(lam-1.0)*0.5 + 1.0 # uncomment to update temp slowly
	vnew = lam*v
	return vnew

def UpdateMP():
	global atomPositions
	global atomVelocities
	#calculate acceleration from forces:
	#print("pre-update positions:")
	if n > 0:
		a = []
		for mol in molecules:
			print ("processing " + mol.name)
			print(atomPositions)
			# TODO calculate values in parallel
			# force fields; truncated Lennard-Jones potential, and EM force via Coulomb's law
			# cumulative method
			#F = -np.array([dLJp(atomPositions,mol,mol.children[atm].indx,sigm[mol.children[atm].tpindex],epsi[mol.children[atm].tpindex],mol.bl[atm]) for atm in range(mol.n)]) # Lennard-Jones
			#F = F-dBEpot(atomPositions,mol,mol.covbonds) # Spring potential. returns array for molecule
			#F = F-dBA(atomPositions,mol,mol.angles) # returns array for molecule
			#F = F-np.array([coul(atomPositions,mol.children[atm].indx,mol.charges[atm],mol,molecules) for atm in range(mol.n)]) # Coulomb
			# rich output:
			lj = -np.array([dLJp(atomPositions,mol,mol.children[atm].indx,sigm[mol.children[atm].tpindex],epsi[mol.children[atm].tpindex],mol.bl[atm]) for atm in range(mol.n)]) # Lennard-Jones
			print("LJ:")
			print(lj)
			ch = -np.array([coul(atomPositions,atm,mol.children[atm].indx,mol.charges[atm],mol,molecules) for atm in range(mol.n)]) # Coulomb
			print("CH:")
			print(ch)
			# bond angles and forces
			bep= -dBEpot(atomPositions,mol,mol.covbonds) # Spring potential. returns array for molecule
			print("bep:")
			print(bep)
			ba = -dBA(atomPositions,mol,mol.angles) # returns array for molecule
			print("ba:")
			print(ba)
			# array of forces on atoms in molecule
			forces= ((lj + bep) + ba) + ch # currently only works for molecules of size 3
			# F=ma     a = F/m     a is extended to include newly calculated accelerations 
			a.extend(np.transpose(np.transpose(force)/mol.mm) for force in forces) #Force->acceleration
			print("acceleration for " + mol.name + ": ")
			print(a)
		
		#print("positions before motion calculations:")
		#print(atomPositions)

		# calculus also allows us to derive velocities and positions from the acceleration
		print("velocities before:")
		print(atomVelocities)
		# verlet integration: atomPositions = 2 * atomPositions - atomPositionsPrev + a
		# requires previous step is saved alongside current step for next step to use. Removes velocity from equations 
		atomVelocities = atomVelocities + np.array(a) * setdt # convert flexible list into numpy array for maths
		print("velocities after:")
		print(atomVelocities)
		atomVelocities = rescaleT(atomVelocities,Temp0) # scale velocities to keep temperature consistent
		atomPositions = atomPositions + atomVelocities * setdt  # this is going to cause problems. I'll have a fixed-rate simulation loop running soon
		#print("updated " + mol.name)
		# split up calculation if dt gets too large to prevent accumulating errors
			
		print("positions at end of UpdateMP():")
		print(atomPositions)

		# boundaries
		if BC == True:
			atomPositions = atomPositions%L
			atomVelocities = 1.0 * atomVelocities
		else:		
			atomPositions, atomVelocities = reflectBC(atomPositions,atomVelocities)

		print("==========END OF AMP UPDATE LOOP==========")
		#total_T_display.text = calculateTemperature(atomVelocities,n)
		#collisions_display.text = str(collisions_count)
		#mouseloc_display.text = str(mouse.collisions)
