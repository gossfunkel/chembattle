from ursina import *
from chempy import Substance
from random import randint
import numpy as np
import atoms as ats

# This will be the game's particle physics engine.
# 				!!! This file currently untested and EXTREMELY janky !!!
# 
# N.B. THIS WILL NOT WORK UNTIL THE UPDATE FUNCTION FROM THE OLD CLASS IS TRANSLATED INTO THIS OR THE ATOMS CLASS
# 	the drawing of the atoms was previously handled by a Molecule object in the AMP class. Atoms will now draw themselves.
# 	
# This physics engine and the networking class will be the main core of the game engine for this project, and if I can 
# 	get a working prototype without wasting the whole summer on it, I think it might be worth keeping up:)

#  ===== === = CLASS DATA = === ===== 

# molecularSim holds the maim simulation data. initialised as a list for loading. Later converted to numpy ndarray
# 	each array on the top level contains the data for a given particle (atom) in the molecular simulation
# 	[0[ typeID ],1[ bonds ],2[ position [x,y,z] ],3[ velocity [x,y,z] ],4[ charge ],5 moleculeID]
molecularSim = []
# stores values for lookup
# 	[0 name ,1[ sigma ],2[ epsilon ],3[ mass ]]
atomTypes 	 = []
# list of models for ursina to draw
drawAtts 	 = []
# list of molecules to calculate bonding  			TODO will this have the bond data in it? probably not
molecules 	 = []

# global values for simulation - can be modified in action
nParticles	= 0 	# number of particles in sim
nMolecules  = 0 	# number of molecules in sim
dimens		= 3 	# number of spacial dimensions
setdt 		= 0.002
simSize		= 15 	# scale of the simulation space in angstroms
simSizeVec	= np.zeros([D])+LL
fixTemp		= 200 	# maintained temperature in K - TODO: CHANGE FOR PRESSURE (walls)

# physical/chemical constants
kR  = 148000.0  #spring potential is usually defined as U = (k/2)(r-r_0)^2. can include /2
kTh = 35300.0  #same explanation as Kr but with bending energy
kb  = 0.8314459920816467 # Boltzmann constant in A^2 D ps^-2 K^-1
kA  = 6.0221409e+26 #Avogadros constant x 1000 (g->kg)
kEq = 1.60217662E-19 #electron charge in coulombs
kQcs = 8.9875517923E9*NA*1E30*ech*ech/1E24 #electrostatic constant in Daltons, electron charges, picosecond, angstrom units

# ===== === = CLASS METHODS = === ===== 

def reflectBC(r,v): 											# update to pass the whole simulation for boundary processing
	newv = 1.0 * v
	newr = 1.0 * r
	for mol in molecules: 										# MOLECULES - CHANGE FOR SIM
		for at in mol.children: 								# ATOMS IN MOLS - CHANGE
			for j in range(dimens): # for each pair of boundaries in a dimension
				if (newr[at.indx,j] < 0): # if beneath 0, reflect
					newr[at.indx,j] = -newr[at.indx,j]
					newv[at.indx,j] = abs(v[at.indx,j])
				if (newr[at.indx,j] > simSizeVec[j]): # if beyond simSize, reflect
					newr[at.indx,j] = 2.0 * simSizeVec[j]-newr[at.indx,j]
					newv[at.indx,j] = -abs(v[at.indx,j])
	return newr,newv 											# TODO update to directly modify positions in simulation array


# Gradient of Lennard-Jones potential on one body (negative force)
def dLJp(molSim,atTypes,atID): # called in loops for every atom in simulation
	#atm = mol.children[atid]
	#pvi = mol.children[atid].indx # position and vector array index of atom

	# if this is the only atom of its type, remove its sigma value to create the temporary array
	#if (atTypes[molSim[atID],3] == 1): 
		#sg=np.delete(sigl,0)
		#eps=
		#print("removed unique type " + tp[mol.tpis[atid]])
	# otherwise index 0 is used for interactions w/ others of its type
	#ep = np.array([epsl])
	# remove epsilon value for self if it is unique
	#if mol.children[atid].isUnique(): ep=np.delete(ep,0) 
	
	#  removal method to keep indices matching commented out
	# create arrays for constants, initially of all type indices
	#sg = [atom[0] for atom in molSim]
	#ep = sg
	#rm = [] # array to remove
	sg = []
	ep = []
	for i in nParticles:
		# populate the array with the values of sigma for each combination of particles
		# n.b. this function is called once for every atom in the sim
		if (molSim[i,5] != molSim[atID,5]): 
			sg.append(atTypes[tpi[i],1,molSim[atID,0]])
			ep.append(atTypes[tpi[i],2, molSim[atID,0]])
			#if (molSim[i,5] == molSim[atID,5]): 
				#rm.append(i)
	#for i in rm:
		#sg = np.delete(sg,i)
		#ep = np.delete(ep,i)
	sixConst = ep*(sg**6)
	twelveConst = 2.0*ep*(sg**12)
	
	# calculate distance vector then absolute distance in each dimension TODO check this is how array indexing works
	drv=molSim[2,axis=1] - molSim[atID,2] # distance vector for every point
	dr=[np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) for a in drv] #absolute distances of that lad

	# truncate to only calculate only for nearby (less than 4 fm away) atoms. May be replaced by partitioning
	if ((dr[i] > 4) for i in dr):
		return 0.0
	else: # calculate dLJP (8-14 from 6-12)
		#print("ep: " + str(ep))
		#print("dr: " + str(dr))
		# TODO
		# 	figure out how the hell to look up the correct sigma and epsilon values for the target interactions
		for i in range(len(sg)):
			r8   = sixConst[i] * (1.0/np.array(dr))**8 # IF THERE'S A DIVIDE BY ZERO, you're calculating particle self-interactions
			r14  = twelveConst[i] * (1.0/np.array(dr))**14
		r8v  = np.transpose(np.transpose(drv)*r8)
		r14v = np.transpose(np.transpose(drv)*r14)
		r8vs = np.sum(r8v,axis=0)
		r14vs= np.sum(r14v,axis=0)
		return 24.0*(r14vs-r8vs)

	#return dLJP
	
# derivative of coulomb potential (negative force)
def coul(molSim, atTypes, atID): # called in loops for every atom in mol, for every mol in molecules 
	#qs = [molSim[i,4] for i in nParticles]
	qs = []
	for i in nParticles:
		if (molSim[i,5] != molSim[atID,5]): 
				qs.append(molSim[i,4])
			# no charge interactions with own molecule
			
	#qs=1.0*np.array([ml.charges for ml in mols])
	qs = 1.0*np.array(qs)
	r = np.delete(molSim[1,axis=1],atID)    # TODO check this is how 2D array indexing works
	drv=r-r[atID] #distance in each dimension
	drv=np.delete(drv,atID,0) #remove ith element (no self-interactions)
	dr=[np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) for a in drv] #absolute distance
	r3=q0*qs*kc*((1.0/np.array(dr))**3.0) # Coulomb's law      TODO Update constants
	FF =np.transpose(np.transpose(drv)*r3) # transpose the distance array, multiply by force, transpose back
	return np.sum(FF,axis=0) # sum the value of axis 0 of the transposed & multiplied array
	#return np.zeros(D) disable charge force

#gradient of bond length potential (negative force)
def dBEpot(molSim): # called once per update
	# 0 initialise a force of 0 on all particles, then fill in those which feel a force
	bps = np.zeros([nParticles,dimens])
	aps = np.zeros([nParticles,dimens])

	# 1 loop through all particle indices in the simulation
	for atID in range(nParticles): 
		#print("calculating bonds for molecule: " + str(mol)
		# position vector of current atom
		rAtom = molSim[atID,2]

		# 2 loop over any bonds in that data entry
		for bond in molSim[atID,1]:
			# load bond length and strength
			# 	bond data structure: [atom,bonded-to,bond-length,gamma(bond-strength)]
			dr0= bond[0] 					# bl
			e0 = molSim[atID,1,1] 			# be

													# TODO DOES THIS NEED TO BE DONE ONE LEVEL UP?
													#	 HOW MANY TIMES ARE BOND ANGLE FORCES BEING CALCULATED?
			aps = dBA(aps,atID,bond[2]) 	# calculate bond angle force around this atom

			# 3 loop over all bonded atoms and calculate forces
			for bondedTo in range(len(bond[3, :])): # iterate through particles bonded to atom
				# print("bonding particle: " + str(atom))
				dr 		   = rAtom - molSim[bnd[bondedTo]] # difference in positions
				#print("dr from dbepot: " + str(dr))
				absdr 	   = np.sqrt(sum(dr**)) 	# squared, summed, then rooted
				dBE 	   = 2.0 * e0 * (adr-dr0) * dr / adr 
				bps[atID] += dBE
	# 4 return the array of forces on all particles in the simulation
	return bps, aps


#gradient of bond angle potential (negative force) 						TODO : REMAKE THIS METHOD FOR NEW DATA STRUCTURE
def dBA(aps,atID,angs): # called once per atom in dBE
	# THIS METHOD IS BUILT FOR WATER-STYLE 3-ATOM MOLECULES 			TODO : ALLOW BONDING TO MORE THAN 2 OTHERS
								# transition metals / coordination complexes will be treated separately, if at all
								# will see if ionic character of bonds pops out on its own or needs some hard-coding
	
	#print("dBA for " + str(atID))
	th00 = angs[0] #equilibrium angle
	e0   = angs[1] #bending modulus
	for i in range(len(angs[2, :])): #loop over particles bonded
		# TODO
		# The original MD code works out angles by going through each particle, and checking an angles list for details
		# 	it then applies a force based on the angle calculated relative to two other particles.
		# 	This reflects sp2-character (trigonal planar) orbital hybridisation geometry, with an implied lone pair 
		# 		on the central atom, and sp hybridisation of the two flanking H atoms.
		# I would like to make this versatile for sp3 (tetrahedral) and sp (linear) geometries, and enable chaining of
		# 	hybridised atoms. Right now, this simply maintains a shape. There is very little chemical logic to it

		# load atomPosition indices from angles array
		a1 = int(angs[i])
		a2 = int(angs[i+1])
		a3 = int(angs[i+2])
		#print("a1: " + str(a1) + ", a2: " + str(a2) + ", a3: " + str(a3))
		# calculate the vectors between the atoms' current positions
		r1 = r[a1] - r[a2] #bond vector 1 (from middle atom to atom 1)
		r2 = r[a3] - r[a2] #bond vector 2 (middle atom to atom 2)
		#print("r values in angle calculation: ra1:" + str(r[a1]) + ", ra2:" + str(r[a2]) + ", ra3:" + str(r[a3]))
		#print("r values in angle calculation: r1:" + str(r1) + ", r2:" + str(r2))
		# calculate absolute lengths of those bonds
		ar1 = abs(np.linalg.norm(r1)) 
		ar2 = abs(np.linalg.norm(r2))
		#print("lengths of bonds in dBA: ar1: " + str(ar1) + ", ar2: " + str(ar2))
		dot = np.dot(r1,r2) #r1 dot r2
		#print("dot product: " + str(dot)) 
		ndot=dot/(ar1*ar2) #normalize dot product by vector lengths i.e. get the cos of angle
		#print("normalised dot product in angle calculation: " + str(ndot))
		th = np.arccos(ndot) #bond angle, theta
		dUdth=-2.0*e0*(th-th00) #-dU/dtheta
		denominator	= np.sqrt(1.0 - ndot**2)
		# if the particle is on the end
		if (i == 2 or i == len(angs)):
			numerator 	= (r2/(ar1*ar2))-(dot/(ar1**3 * ar2 * 2.0))
			dUdr 		= dUdth * numerator/denominator
			aps[i+1] 	+= dUdr
		# else it is in the middle
		else:
			#print (" denominator in angle calculation: " + str(denominator))
			numerator	= (dot * r1 / ar1**2 + dot * r2 / ar2**2 - (r2 + r1))/(ar1*ar2)
			dUdr		= dUdth * numerator/denominator
			aps[i+1]	+= dUdr
			aps[i]		+= dUdr
			aps[i-1]	+= dUdr
	return aps

def CreateMol(molecularSim, atomTypes, molName, numMolsCreate, player, location):

	atfact  = ats.AtomFactory() # factory method for generating atoms
	newatts = [] 			    # instantiate list for new atoms in molecule

	# interpret molName
	newmol 		 = Substance.from_formula(molName) # get a new chempy Substance object from the name
	comp 		 = newmol.composition 				 # get the elemental composition as a dictionary
	bondLengs 	 = [] 																		# TODO find bond lengths
	bondEnergies = [] 																		# TODO find bond energies
	bondAngles	 = [] 																		# TODO find bond angles
	for elementindex,quant in comp.items(): # iterate through the members of the composition dictionary
		for i in range(quant): # loop for quantity
			# generate atoms of given type
			newt = atfact.createAtom(elementindex,player)
			newatts.append(newt) # add new atom to list
			#print("created " + newt.name)
	
	# initialise array for simulation particle data entries for the new atoms in the molecule
	newDataArray = np.empty(len(newatts))
	molSize 	 = len(newatts)
	sumcharge 	 = 0
	elnegal = []
	sumelecneg = 0

	# load each atom
	for newAtom in range(molSize):
		nAtomTypes  = len(atomTypes)
		# initialise new index as the length of the atomTypes array (1 above highest valid index)
		# 	moleculeSim atom data format: [typeIndex,bonds,angles,position,velocity,charge]
		newAtomData = [nAtomTypes,[],[],[],[],0]
		newat 		= newatts[newAtom]
		# check types: search typearray for atom name
		for typ in range(nAtomTypes):
			# if the name in the given type array matches, save its index to the new atom's data array
			if (atomTypes[typ,0] == newat.name):
				newAtomData[0] = typ
		
		if (newAtomData[0] == nAtomTypes): # create new type and add to atomTypes 
			# 1 save type array index for data
			newAtomData[0] = tpindex
			# 2 load values - name, sigma(i), epsilon(i),mass,population=1 for newly initialised type
				# atomType format:	[name,[ sigma ],[ epsilon ],mass]
			atomTypes.append(newat.name,[newat.sig],[newat.eps],newat.mass],1)
		else: 
			# type already loaded, index value already copied
			# increment population of atomType
			atomTypes[tpindex,4] = atomTypes[tpindex,4] + 1
		# load bond lengths, energies, and angles
			#covbonds structure: [indx,bonded-to,bond-length,gamma(bond-strength)]
		if (molSize > 1):
			newAtomData[newAtom,1] = [bondLengs[newAtom],kR]]
			if newAtom < (molSize - 1):
																				# TODO program in variable bonding configurations
				newAtomData[2] = [newAtom-1,newAtom,newAtom+1,bondAngles[newAtom],kTh]
		# load initial position values
		if (player == 0 and location == None): newAtomData[3] = Vec3(-30+randint(-5,5),-12+randint(-1,1),-20+randint(-5,5))
		elif (player == 1 and location == None): newAtomData[3] = Vec3(30+randint(-5,5),-12+randint(-1,1),20+randint(-5,5))
		elif (player == 0 and location != None): newAtomData[3] = Vec3(-30+randint(-5,5),-12+randint(-1,1),-20+randint(-5,5)) + location
		elif (player == 1 and location != None): newAtomData[3] = Vec3(30+randint(-5,5),-12+randint(-1,1),20+randint(-5,5)) + location

		if (newAtom > len(bondLengs)-1):
			newAtomData[newAtom,1] = bondLengs[newAtom-1]
			if (newAtom > 0):
				# increment temporary position variable by bond length to give smoother load-in
				newAtomData[3] = np.array(newAtomData[3]) + np.array([1.0/bondLengs[newAtom]**3 , 1.0/bondLengs[newAtom]**3 , 0])
					# in future, this will be replaced with / followed by Maxwell-Boltzmann to match initial velocities to Temp
		
		# tell atom game objects where they will be rendered
		newat.world_location = Vec3(newAtomData[3])

		# initialise velocity to 0
		newAtomData[4] = [np.zeros(dimens)]

		# do charge calculations
		formalCharges.append(newat.charge)
		sumcharge += newat.charge
		elnegal.append(newat.electroNegAllen)
		sumelecneg += newat.electroNegAllen # ========== END OF FOR LOOP - THRU ATOMS ============
	# calculate charge distribution
	# hacky temporary assigning of partial charge based on electronegativity
	# calculate partial charges using Allen electronegativity and formal charge contributions
	formalCharges = np.array(formalCharges)
	# might use this? if i don't ditch all of this for a better method
	#sortedChrg = np.array(chrg)
	#sortedChrg.sort()
	#diffChrg = sortedChrg[-1] - sortedChrg[0]

	# find the difference between the largest and smallest electronegativities in the atom
	sortedElnegal = np.array(elnegal)
	sortedElnegal.sort()
	diffNeg = sortedElnegal[-1] - sortedElnegal[0]
	# go through atoms and calculate partial charge as a function of electronegativity
	for atomnum in range(molSize):
		if (formalCharges[atomnum] < 0): 	# assign sign based on charge of free atom
			diffNeg *= -1
		if sumcharge == 0:  # this molecule has a formal charge of 0
			newAtomData[atomnum,5] = diffNeg * (elnegal[i] / sumelecneg) # latter term represents ratio of electronegativity in molecule
		else: # this molecule has an overall charge on it; it's an ion
			newAtomData[atomnum,5] = sumcharge + diffNeg * (elnegal[i] / sumelecneg) # latter term represents ratio of electronegativity in molecule
			
		# AMBER uses e = sum(qiqj/rij) so could maybe rearrange. distance could define how much electrons can split charge
		# where at a distance beyond the maximum bonding distance the charge on an atom reaches the formal charge of a free atom
		# this might have to be saved to see if I can figure out a holistic approach to electrons

	# load new data into simulation 
	molSim = np.ndarray(listSim)
	nMolecules 	+= 1
	listSim 	 = molecularSim.toList()
	listSim.append([newData for newData in newAtomData]) # TODO add some exception raising here if something goes wrong
	# pass models to physics class
	drawAtts.append([newat for newat in newatts])
	return molSim # not sure if this makes sense, dont know how to call this and not have scope issues
							# since molecularSim isn't being used as a global variable, I should make a setter to call this
							# complex routine from outside the simulation, to give more moderation on who can call it

# decision logic for when to construct a molecule at code request
# 	Just because CreateMol can be called, doesn't mean it should execute. This function controls that and
# 		helps manage scope.
def RequestCreateMol(atomTypes, molName, numMolsCreate, player=0, location=None):
	global molecularSim
										# TODO move things like positioning of atoms etc to this method
	molecularSim = CreateMol(atomTypes, molName, numMolsCreate, player, location)
	# throw some exceptions if something goes wrong loading the molecules; expect a True return
	return True

# ===== === = UPDATE FUNCTION = === =====

def Update():
	# reinitialise array of forces	==========		TODO this is sooooo cycle inneficient, pls redo
	forces = np.zeros(nParticles,dimens)
	for atom in range(nParticles):
		# currently loading to variables rather than accumulating F for debugging purposes
		lj = -np.array([dLJp(molecularSim, atomTypes, atom)]) # Lennard-Jones
		#print("LJ: + str(lj))
		ch = -np.array([coul(molecularSim, atomTypes, atom)]) # Coulomb
		#print("CH:" + str(ch))
		forces[atom] = [lj,ch]

	# 							==========
	# bond angles and forces
	bep= -dBEpot(molecularSim) # Spring potential. returns array for molecule
	#print("bep: + str(bep))
	ba = -dBA(molecularSim) # returns array for molecule
	#print("ba:" + str(ba))
	# array of forces on atoms in molecule
	forces= [((lj[i] + bep[i]) + ba[i]) + ch[i] for i in range(nParticles)]
	# F=ma     a = F/m     			a is extended to add all newly calculated accelerations 
	a.extend(np.transpose(np.transpose(forces)/molSim[MASSES])) #Force->acceleration 		TODO store a mass value for each particle
	#print("acceleration at end of calculations: " + str(a))
	# 							==========

	# calculus also allows us to derive velocities and positions from the acceleration
	# CONSIDER: verlet integration: r[t+1] = 2 * r[t] - r[t-1] + a
	# 	requires previous step is saved alongside current step for next step to use. Removes velocity from equations 
	# 	would take the place of velocity array in the simulation, could make for more reliable networking/conflict resolution
	molecularSim[4] = molecularSim[4] + np.array(a) * setdt # convert flexible list of accellerations into numpy array for maths
	molecularSim[4] = rescaleT(molecularSim[4],Temp0) # scale velocities to keep temperature consistent TODO CHANGE FOR PRESSURE FROM WALLS

	# UPDATE PARTICLE POSITIONS ==========
	molecularSim[3] = molecularSim[3] + molecularSim[4] * setdt
	# still need to figure out frame-rate consistency. setdt is required to get reliable behaviour/physics

	# process boundary condition: reflect off of walls. 		TODO maintain constant pressure so that temp can vary
	molecularSim[3] , molecularSim[4] = reflectBC(molecularSim[3] , molecularSim[4])
		
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
