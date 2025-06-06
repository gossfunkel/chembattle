from ursina import *
from chempy import Substance

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
# This physics engine and the networking class will be the main core of the game engine for this project, and if I can 
# get a working prototype without wasting the whole summer on it, I think it might be worth keeping up:)
# '''

# holds the simulation data
# initialised as a list for loading. Later converted to numpy ndarray
# each array on the top level contains the data for a given particle (atom) in the molecular simulation
# [0[ typeID ],1[ bonds ],2[ position [x,y,z] ],3[ velocity [x,y,z] ],4[ charge ],5 moleculeID]
molecularSim = []
# stores values for lookup
# [0 name ,1[ sigma ],2[ epsilon ],3[ mass ]]
atomTypes 	 = []
# list of models to draw
drawAtts 	 = []

# global values for simulation - can be modified in action
nParticles	= 0 	# number of particles in sim
nMolecules = 0
dimens		= 3 	# number of spacial dimensions
setdt 		= 0.002
simSize		= 15 	# scale of the simulation space in angstroms
simSizeVec	= np.zeros([D])+LL
fixTemp		= 200 	# maintained temperature in K - TODO: CHANGE FOR PRESSURE (walls)

kR  = 148000.0  #spring potential is usually defined as U = (k/2)(r-r_0)^2. can include /2
kTh = 35300.0  #same explanation as Kr but with bending energy
kb  = 0.8314459920816467 # Boltzmann constant in A^2 D ps^-2 K^-1
kA  = 6.0221409e+26 #Avogadros constant x 1000 (g->kg)
kEq = 1.60217662E-19 #electron charge in coulombs
kQcs = 8.9875517923E9*NA*1E30*ech*ech/1E24 #electrostatic constant in Daltons, electron charges, picosecond, angstrom units

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

def CreateMol(molecularSim, atomTypes, molName, numMolsCreate, player=0, location=None):

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
	listSim = molecularSim.toList()
	listSim.append([newData for newData in newAtomData])
	molecularSim = np.ndarray(listSim)
	# pass models to physics class
	drawAtts.append([newat for newat in newatts])
	return molecularSim # not sure if this makes sense

def Update():
	for atom in nParticles:
		lj = -np.array([dLJp(molecularSim, atomTypes, atom)]) # Lennard-Jones
		#print("LJ:")
		#print(lj)
		ch = -np.array([coul(molecularSim, atomTypes, atom)]) # Coulomb
		#print("CH:")
		#print(ch)
		# bond angles and forces
		bep= -dBEpot(molecularSim) # Spring potential. returns array for molecule
		#print("bep:")
		#print(bep)
		ba = -dBA(molecularSim) # returns array for molecule
		#print("ba:")
		#print(ba)
		# array of forces on atoms in molecule
		F= ((lj + bep) + ba) + ch # currently only works for molecules of size 3
		# F=ma     a = F/m     a is extended to include newly calculated accelerations 
		a.extend(np.transpose(np.transpose(F)/mol.mm)) #Force->acceleration
		#print("acceleration for " + mol.name + ": ")
		#print(a)
	# calculus also allows us to derive velocities and positions from the acceleration
	#print("velocities before:")
	#print(atomVelocities)
	# verlet integration: atomPositions = 2 * atomPositions - atomPositionsPrev + a
	# requires previous step is saved alongside current step for next step to use. Removes velocity from equations 
	molecularSim[4] = molecularSim[4] + np.array(a) * time.dt # convert flexible list into numpy array for maths
	molecularSim[4] = rescaleT(molecularSim[4],Temp0) # scale velocities to keep temperature consistent
	molecularSim[3] = molecularSim[3] + molecularSim[4] * time.dt  # this is going to cause problems. I'll have a fixed-rate simulation loop running soon

	molecularSim[3] , molecularSim[4] = reflectBC(molecularSim[3] , molecularSim[4])
		