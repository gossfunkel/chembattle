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
	# 	[0[ typeID ],1[ bonds&angles ],2[ position [x,y,z] ],3[ velocity [x,y,z] ],4[ charge ],5 moleculeID,6 mass]
molecularSim = []
# stores values for lookup
# 	[0 name ,1[ sigma ],2[ epsilon ]]
atomTypes 	 = []
# list of models for ursina to draw
drawAtts 	 = []
# list of molecules to use for calculating bond potentials
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
kB  = 0.8314459920816467 # Boltzmann constant in A^2 D ps^-2 K^-1
kA  = 6.0221409e+26 #Avogadros constant x 1000 (g->kg)
kEq = 1.60217662E-19 #electron charge in coulombs
kQcs = 8.9875517923E9*NA*1E30*ech*ech/1E24 #electrostatic constant in Daltons, electron charges, picosecond, angstrom units

# ===== === = CLASS METHODS = === ===== 

def reflectBC(r,v):
	newv = 1.0 * v
	newr = 1.0 * r
	for atom in range(nParticles):
		for j in range(dimens): # for each pair of boundaries in a dimension
			if (newr[atom,j] < 0): # if beneath 0, reflect
				newr[atom,j] = -newr[atom,j]
				newv[atom,j] = abs(v[atom,j])
			if (newr[atom,j] > simSizeVec[j]): # if beyond simSize, reflect
				newr[atom,j] = 2.0 * simSizeVec[j]-newr[atom,j]
				newv[atom,j] = -abs(v[atom,j])
	return newr,newv

def rescaleT(v):
	# mass / v^2
	kineticE   = 0.5 * sum(sum(molSim[:,6] * np.transpose(v**)))
	avKineticE = KineticE / nParticles
	tNow 	   = (2.0 / 3) * avKineticE / kB
	lam 	   = np.sqrt(temp0 / tNow)
	#lam = (lam - 1.0) * 0.5 + 1.0 # uncomment to update temp slowly
	return (lam * v) # adjusted velocity

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

# Gradient of Lennard-Jones potential on one body (negative force)
def dLJp(atID): # called in loops for every atom in simulation
	# convenience variable
	atTpID = molecularSim[atID,0]
	#pvi = mol.children[atid].indx # position and vector array index of atom

	# create arrays for constants for the type of the atom being calculated for
	sg = [atomTypes[:,1,atTpID].copy()]
	ep = [atomTypes[:,2,atTpID].copy()]
	
	# calculate distance vector then absolute distance in each dimension 
	drv  = molecularSim[:, 2].copy() - molecularSim[atID,2].copy() 	# distance vector for every point
	dr   = [np.sqrt((a[0]**)+(a[1]**)+(a[2]**)) for a in drv] 		# absolute distances of those vectors
	#print("ep: " + str(ep) + " || dr: " + str(dr))
	dLJP = np.zeros(dimens) 									 	# force vector on atom being calculated for

	# truncate to calculate only for nearby atoms (less than 3 angstroms away). May be replaced by partitioning
	if ((dr[i] < 3) for i in dr): 				# calculate dLJP (8-14 from 6-12) force from distance array member
		# NOTE: IF THERE'S A DIVIDE BY ZERO, you're calculating particle self-interactions
		r8   = ep[i]*(sg[i]**6) * (1.0/np.array(dr))**8 
		r14  = 2.0*ep[i]*(sg[i]**12) * (1.0/np.array(dr))**14
		r8v  = np.transpose(np.transpose(drv)*r8)
		r14v = np.transpose(np.transpose(drv)*r14)
		r8vs = np.sum(r8v,axis=0)
		r14vs= np.sum(r14v,axis=0)
		dLJP += 24.0*(r14vs-r8vs)
	# return the total force vector on the atom
	return dLJP
	
# derivative of coulomb potential (negative force)
def coul(atID): # called once per particle, calculates field effect as a whole
	# TBC: Will see if ionic character of bonds pops out or needs some hard-coding
	#qs = [molecularSim[i,4] for i in nParticles]
	qs = [molecularSim[i,4].copy() for i in nParticles if (molecularSim[i,5] != molecularSim[atID,5])]
			# no charge interactions with own molecule
	qs 			= 1.0 * np.array(qs)
	r 			= np.delete(molecularSim[:, 1].copy(),atID)
	drv 		= r - r[atID] #distance in each dimension
	drv 		= np.delete(drv,atID,0) #remove ith element (no self-interactions)
	dr 			= [np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) for a in drv] #absolute distance
	r3 			= molecularSim[atID,4].copy() * qs * kQcs * ((1.0 / np.array(dr)**3.0) # Coulomb's law      TODO Update constants
	forcefield  = np.transpose(np.transpose(drv) * r3) # transpose the distance array, multiply by force, transpose back
	#return np.zeros(dimens) # uncomment to disable EM field force
	return np.sum(forcefield,axis=0) # sum the value of axis 0 of the transposed & multiplied array

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
			
			# 3 loop over all bonded atoms and calculate forces
			for bondedTo in range(len(bond[3, :])): # iterate through particles bonded to atom
				#print("bonding particle: " + str(atom))
				dr 		   = rAtom - molSim[bnd[bondedTo]] # difference in positions
				#print("dr from dbepot: " + str(dr))
				bps[atID] += 2.0 * e0 * (np.sqrt(sum(dr**)) - dr0) * dr / adr 
				# 4 Calculate bond angle forces (once per bonded neighbour)
				aps = dBA(aps,atID,bnd[bondedTo],bond[2,0],bond[2,1])
				# because this calculates through the neighbours, it is not suited to calculating the force
				# on atoms in-between other atoms. I may write a new procedure for this; dBA is getting a bit mad

			#print("calculating dBA for " + str(atID))
			#aps = dBA(aps,atID,bond[2]) 	# calculate bond angle force around this atom
			)
	# 5 return the array of forces on all particles in the simulation
	return bps, aps

#gradient of bond angle potential (negative force) 						TODO : REMAKE THIS METHOD FOR NEW DATA STRUCTURE
def dBA(aps,atomID,neighbour,th00,e0): # called once per bonded atom to each atom in dBE
	# this calculates bond angle forces for every pair of bonds	around a shared/middle/mediating/axial atom
		# So if a particle has multiple bonds, forces will need to be calculated multiple times.
		# If angles don't properly match, geometry may get fucky.
		# It updates an array of force on atoms from bond angles (aps)
	# The original MD code works out angles by going through each particle, and checking an angles list for details
		# 	it then applies a force based on the angle calculated relative to two other particles.
		# 	This reflects sp2-character (trigonal planar) orbital hybridisation geometry, with an implied lone pair 
		# 		on the central atom, and sp hybridisation of the two flanking H atoms.
		# I would like to make this versatile for sp3 (tetrahedral) and sp (linear) geometries, and enable chaining of
		# 	hybridised atoms. Right now, this simply maintains a shape. There is very little chemical logic to it. This 
		# 	could continue, but I want to include lone pairs- for now, sp3 hybrids will simply work by stacking bonds on an
		# 	atom, and sp hybrids will break the routine (ndot of opposing vectors is -1, results in dividing by 0)
	# The new methodology of this method is that it is called once per bonded atom. Angle forces are exerted by neighbours
		# of each neighbour, so this loop iterates through all neighbours of the neighbour and calculates forces to maintain
		# the bond angles thereof. It outputs only the force on each side-particle; I have not yet implemented a solution to
		# the lack of force on the atoms acting as an intermediary/middle/axis between other bonded neighbours. The 
		# numerator in the differentiation for the force on the middle particle is given as: 
			#	numerator		  = (dot * r1 / ar1**2 + dot * r2 / ar2**2 - (r2 + r1)) / (ar1*ar2)
	# N.B.: due to the choice to calculate for every particle, calculations will need to be done perhaps more than they need 
		# to. Given particles i, j, and k, where i is bonded to j and j is bonded to k, this will calculate the bond angle ijk
		# twice; once when calculating for i, and once when calculating for k. In future, I'll probably have to rewrite this
		# to be more efficient with the computations. I haven't yet figured out how to do this.
											
											# 	TODO : MAKE ROUTINE CHECK IF BOND ANGLE HAS ALREADY BEEN CALCULATED BEFORE
									# 	this will be necessary to prevent each atom's bond array generating duplicate forces

	# 1 set values for angle to enforce
	#th00 = angs[0] #equilibrium angle
	#e0   = angs[1] #bending modulus
	if (th00 > 170 and th00 < 190): # sp character - linear
		# code for calculating bond angle force for linear bonds
		aps[atID] = 0.0
	elif(th00 > 110 and th00 < 130): # sp2 character - trigonal planar
		# 2 loop over particles bonded to neighbour -1 (no two-member angles! angles need 3+ points)
		for angledAtom in range(len(moleculeSim[neighbour,1])) if (moleculeSim[neighbour,1,angledAtom] != atID): 
			# 2.1 load indices from angles array
			a3 = int(moleculeSim[neighbour,1,angledAtom]])
			#print("a1: " + str(a1) + ", a2: " + str(a2) + ", a3: " + str(a3))

			# 2.2 calculate the vectors between the atoms' current positions
			r1 = r[atomID] - r[neighbour] #bond vector 1 (from middle atom to atom currently being calculated for)
			r2 = r[a3] - r[neighbour] #bond vector 2 (middle atom to other atom)
			#print("r values in angle calculation: ra1:" + str(r[a1]) + ", ra2:" + str(r[a2]) + ", ra3:" + str(r[a3]))
			#print("r values in angle calculation: r1:" + str(r1) + ", r2:" + str(r2))

			# 2.3 calculate absolute lengths of those vectors, their dot product, then normalise and arccos to get angle theta
			ar1 = abs(np.linalg.norm(r1)) 
			ar2 = abs(np.linalg.norm(r2))
			#print("lengths of bonds in dBA: ar1: " + str(ar1) + ", ar2: " + str(ar2))
			dot = np.dot(r1,r2) 						# dot product of distance vectors between atoms in molecule
			#print("dot product: " + str(dot)) 
			ndot=dot/(ar1*ar2) 							# normalize dot product by vector lengths i.e. get the cos of angle
			#print("normalised dot product in angle calculation: " + str(ndot))
			th = np.arccos(ndot) 						# bond angle, theta
			#if (th == 180): raise Exception("180 degree bond angles not supported by dBA routine yet!")
			
			# 2.4 check atom position in angle 							TODO maths stuff to get rid of these if statements
			#if (i == 2 or i == len(angs)):  			# if the particle is on the end of the molecule 
														# calculation for differentiation - linear end atoms
			numerator 		  = (r2/(ar1*ar2))-(dot/(ar1**3 * ar2 * 2.0))
			#else: 										# else it is in the middle
			# TODO can I just remove this? I guess the mediating atom also feels a force from the angle that I should include
			#											# calculation for differentiation - axial atoms
			#	numerator		  = (dot * r1 / ar1**2 + dot * r2 / ar2**2 - (r2 + r1)) / (ar1*ar2)
			# 2.5 differentiate to calculate force
			dUdth = -2.0 * e0 * (th-th00) # - dU/dtheta
			denominator	= np.sqrt(1.0 - ndot**2)
			#print (" denominator in angle calculation: " + str(denominator))

			# 2.5.1 force found as derivative of energy over distance (force field) i.e. dUdr
			aps[atomID] += dUdth * numerator/denominator
	return aps

def CreateMol(molName, numMolsCreate, player, position): 							# BROKEN : BONDS, ANGLES
	# 1 initialise variables
	atfact  = ats.AtomFactory() # factory method for generating atoms
	newatts = [] 			    # instantiate list for new atoms in molecule
	atTypes = atomTypes
	molID   = nMolecules + 1
	
	# 2 interpret molName
	newmol 		 = Substance.from_formula(molName) # get a new chempy Substance object from the name
	comp 		 = newmol.composition 				 # get the elemental composition as a dictionary
	bondLengs 	 = [] 																		# TODO find/assign bond lengths
	bondEnergies = [] 																		# TODO find/assign bond energies
	bondAngles	 = [] 																		# TODO find/assigm bond angles
	for elementindex,quant in comp.items(): # iterate through the members of the composition dictionary
		for i in range(quant): # loop for quantity
			# generate atoms of given type and add new atom to list
			newatts.append(atfact.createAtom(elementindex,player)) 
	
	# 3 initialise array for simulation particle data entries for the new atoms in the molecule
	newDataArray = np.empty(len(newatts)).toList()
	molSize 	 = len(newatts)
	sumcharge 	 = 0
	elnegal = []
	sumelecneg = 0

	# 4 load each atom
	for newAtom in range(molSize):
		nAtomTypes  = len(atTypes)
		# 4.1 initialise new index as the length of the atomTypes array (1 above highest valid index)
		# 			and load the mass of the particle													TODO MASS !!!!!!!!!!!
		# 	moleculeSim atom data format: [typeIndex,bonds,position,velocity,charge,moleculeID,mass]
		newAtomData = [nAtomTypes,[],[],[],[],molID,MASS]
		newat 		= newatts[newAtom]
		# 4.2 check types: search typearray for atom name
		# if the name in the given type array matches, save its index to the new atom's data array
		newAtomData[0] = [typ for typ in range(nAtomTypes) if (atTypes[typ,0] == newat.name)]
		
		if (newAtomData[0] == nAtomTypes): # create new type and add to atomTypes 
			# 4.2.1 save type array index for data
			newAtomData[0] = tpindex
			# 4.2.2 load values - name, sigma(i), epsilon(i),population=1 for newly initialised type
				# atomType format:	[name,[ sigma ],[ epsilon ],mass]
			atTypes.append(newat.name,[newat.sig],[newat.eps]],1)
		else: 
			# type already loaded, index value already copied
			# increment population of atomType
			atTypes[tpindex,4] = atTypes[tpindex,4] + 1
		# 4.3 load bond lengths, energies, and angles
			#covbonds structure: [indx,bonded-to,bond-length,gamma(bond-strength)]
		if (molSize > 1):
			#newAtomData[1] = [bondLengs[newAtom],kR]]
			if newAtom < (molSize - 1):
																			# TODO HOW TO TELL NEW ATOMS HOW TO BOND?
				for bonds in atomBonds:
					# TODO bond length and atoms to bond is information that this constructor needs
					# 	this information will have to be loaded with the molecule and assigned here
					newAtomData[1].append([bondLength,kTh,newAtom-1,newAtom,newAtom+1,bondAngles[newAtom],])
		# 4.4 load initial position values
		newAtomData[2] = position
		# 4.5 tell atom game objects where they will be rendered
		newat.world_location = Vec3(newAtomData[2])
		# 4.6 increment position for next atom generated 					TODO !!! NEED NEW BOND LENGTHS VALUES/ARRAY
		position = np.array(newAtomData[2]) + np.array([1.0/bondLengs[newAtom]**3 , 1.0/bondLengs[newAtom]**3 , 0])
		if (i%2==0):
			#self.children[i].rotation_y = 180
			position[1] += 1
		else:
			position[2] -= 1
		# temporarily extend bond info along whole molecule 	TODO this seems wrong
		#if (newAtom > len(bondLengs)-1):
		#	newAtomData[newAtom,1] = bondLengs[newAtom-1]

		# 4.7 initialise velocity to 0
		newAtomData[3] = [np.zeros(dimens)]

		# 4.8 do charge calculations
		formalCharges.append(newat.charge)
		sumcharge += newat.charge
		elnegal.append(newat.electroNegAllen)
		sumelecneg += newat.electroNegAllen 

		newDataArray.append([newAtomData])# ========== END OF FOR LOOP THRU NEW ATOMS ============
	
	# 5 calculate charge distribution
								# hacky temporary assigning of partial charge; calculating partial charges
								# 	using Allen electronegativity and formal charge contributions
								#	 might use this? if i don't ditch it for a better method
	formalCharges = np.array(formalCharges)
	# 5.1 find the difference between the largest and smallest electronegativities in the atom
	sortedElnegal = np.array(elnegal)
	sortedElnegal.sort()
	diffNeg = sortedElnegal[-1] - sortedElnegal[0]
	# 5.2 go through atoms and calculate partial charge as a function of electronegativity
	for atomnum in range(molSize):
		# assign sign based on charge of free atom
		diffNeg = formalCharges[atomnum]/formalCharges[atomnum]
		newDataArray[atomnum,5] = sumcharge + diffNeg * (elnegal[atomnum] / sumelecneg) 
		# latter term represents ratio of electronegativity in molecule. 
			
		# AMBER uses e = sum(qiqj/rij) so could maybe rearrange that. distance could define how much electrons can split charge
		# where at a distance beyond the maximum bonding distance the charge on an atom reaches the formal charge of a free atom
		# this might have to be saved to see if I can figure out a holistic approach to electrons

	# 6 load new data into simulation 
	nMolecules 	+= 1
	listSim 	 = molecularSim.toList()
	listSim.extend(newDataArray) 			# TODO add some exception raising here if something goes wrong
	# 7 pass models to list to render
	drawAtts.extend(newatts)
	return np.ndarray(listSim), atTypes # not sure if this makes sense, dont know how to call this and not have scope issues
							# since molecularSim isn't being used as a global variable, I should make a setter to call this
							# complex routine from outside the simulation, to give more moderation on who can call it

# decision logic for when to construct a molecule at code request. Should return True if succesful.
		# 	Just because CreateMol can be called, doesn't mean it should execute. This function controls that and
		# 		helps manage scope, and adds some game logic packaging
def RequestCreateMol(molName, numMolsCreate, player=0, location=None):
	# 	TODO allow passing multiple types of molecule at once with *variable passing

	# 1 calculate molecule spawn location based on which player is spawning the molecules
	if (location == None): newAtomData[3] = Vec3((player+1)*18 + randint(-5,5),
															 5 + randint(-1,1),
							 (100*(player+1) - 100 + 8*player) + randint(-5,5))
	else: newAtomData[3] = Vec3((player+1)*18 + randint(-5,5),
											5 + randint(-1,1),
			(100*(player+1) - 100 + 8*player) + randint(-5,5)) + location

	# 2 create the new molecule and update the simulation to include it
	global molecularSim
	global atomTypes
	molecularSim, atomTypes = CreateMol(atomTypes, molName, numMolsCreate, player, location)
	# 3 throw some exceptions if something goes wrong loading the molecules; expect a True return to this function
	if type(molecularSim) == None: 
		molecularSim = []
		raise Exception("MOLECULAR SIM NOT LOADED PROPERLY; RESETTING")

	return True

# ===== === = UPDATE FUNCTION = === =====

def Update():
	# 0 reinitialise array of forces	==========		TODO this is sooooo cycle inneficient, pls redo
	forces = np.zeros(nParticles,dimens)
	# 1 iterate through particles in simulation
	forces[atom] = [-np.array([dLJp(atom)]),-np.array([coul(atom)]) for atom in range(nParticles)]

	# 							==========
	# 2 calculate forces maintaining bond angles and lengths
	bep= -dBEpot(molecularSim) # Spring potential. returns array for molecule
	#print("bep: + str(bep))
	ba = -dBA(molecularSim) # returns array for molecule
	#print("ba:" + str(ba))
	# 3 calculate array of forces on all atoms - LJ -> bonds -> EM field
	forces= [((forces[i,0] + bep[i]) + ba[i]) + forces[i,1] for i in range(nParticles)]
	# 4 F=ma     a = F/m     			a is extended to add all newly calculated accelerations 
	a = np.transpose(np.transpose(forces)/molSim[:,6]) #Force->acceleration 
	#print("acceleration at end of calculations: " + str(a))
	# 							==========

	# 5 calculus also allows us to derive velocities and positions from the acceleration
	# Euler method. This is unstable and creates increasing errors, especially for anything other than the shortest timesteps
		# CONSIDER: verlet integration: r[t+1] = 2 * r[t] - r[t-1] + a
		# 	requires previous step is saved alongside current step for next step to use. Removes velocity from equations 
		# 	would take the place of velocity array in the simulation, could make for more reliable networking/conflict resolution
		# CONSIDER: rk4 (runge-kutte)
		# CONSIDER: Beeman predictor-corrector
		# CONSIDER: Conditional methodology
	molecularSim[:,3] = molecularSim[:,3] + np.array(a) * setdt # convert flexible list of accellerations into numpy array for maths
	molecularSim[:,3] = rescaleT(molecularSim[:,3]) # scale velocities to keep temperature consistent TODO CHANGE FOR PRESSURE FROM WALLS

	# 6 UPDATE PARTICLE POSITIONS ==========
	molecularSim[:,2] = molecularSim[:,2] + molecularSim[:,3] * setdt
	# still need to figure out frame-rate consistency. setdt is required to get reliable behaviour/physics

	# 7 process boundary condition: reflect off of walls. 		TODO maintain constant pressure so that temp can vary
	molecularSim[:,2] , molecularSim[:,3] = reflectBC(molecularSim[:,2], molecularSim[:,3])
		
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
