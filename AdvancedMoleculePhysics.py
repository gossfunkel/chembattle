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
# ''' TODO : REWRITE ALL THE FUNCTIONS BELOW TO FIT THE NEW STRUCTURE OF THE ARRAYS
#
# This is the game's physics engine.
# 
# I am in the process of rewriting a basic Lennard-Jones, Spring Bonding and Electric molecular dynamics sim. I've integrated Ursina, but not yet
# finished creating the core atom loading and processing algorithms.
# 
# In short, don't expect this code to run reliably. The rest of the game isn't built up to keep up with these changes yet; hence 
# splitting this off into a new physics file to keep the original sandbox going in parallel. This physics engine and the multiplayer
#  are the main core of the game engine for this project, and if I can get a working prototype without wasting the whole summer on 
# it, I think it might be worth keeping up. 
# '''

# global values for simulation - can be modified in action
n 			= 0 	# number of particles in sim
D 			= 3 # number of spacial dimensions
setdt 		= 0.002
LL 			= 20 	# scale of the simulation space
L 			= np.zeros([D])+LL
Temp0 		= 200 	# maintained temperature - TODO: CHANGE FOR PRESSURE (walls)
#dragstate 	= False # bool value to pause the sim while dragging a molecule around

# constants
Kr  = 148000.0  #spring potential is usually defined as U = (k/2)(r-r_0)^2. can include /2
Kth = 35300.0  #same explanation as Kr but with bending energy
kb  = 0.8314459920816467 # Boltzmann
NA  = 6.0221409e+26 #Avogadros constant x 1000 (g->kg)
ech = 1.60217662E-19 #electron charge in coulombs
kc  = 8.9875517923E9*NA*1E30*ech*ech/1E24 #electrostatic constant in Daltons, electron charges, picosecond, angstrom units

# arrays of generated particles
# TODO: Make 3D so that each molecule can have its own array of positions nested, to lookup by mol and child indices
atomPositions 	= [] 	# will be later initialised into a fast numpy array of vectors
atomVelocities 	= [] 	# will be later initialised into a fast numpy array of vectors
molecules 		= []	# list storing molecule objects
tp 				= [] 	# list tracking species of atom (element/isotope)
masses 			= [] 	# list tracking masses
a 				= [] 	# array for use by update function to index arrays of molecules' components' accelerations
sigm 			= []	# array for population with sigma (balanced potential) LJP values for loaded atoms
epsi 			= []	# array for population with epsilon (potential well) LJP values for loaded atoms
#chrgs 			= [] 	# array for use by coul function to index all charges

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
	def __init__(self, nam, position, indx, chldatoms, tpis, velocity=np.zeros(D), bl=[1], th0=109.47*np.pi/180.0):
		if (len(chldatoms) == 0): return None # don't construct an empty molecule; what is that??? how would it even work?!
		self.n 			  = len(chldatoms)
		super().__init__(world_position=position)
		self.nam 		  = nam
		self.visible_self = False # molecule isn't visible; atoms are visible
		#if self.n > 0:
		self.children = chldatoms # the children of the molecule object are the atoms that constitute the molecule
		#print(self.children)
		self.indx 		  = indx # the index of the first child atom in the position and velocity arrays
		# self.bl 		  = *bl 	# 0.9611 for water -- bond lengths stored in the arrays stored in covbonds
		self.mm 		  = [] 	# array storing masses
		
		# construction lists
		covbonds		  = {}	# dict tracking covalent bonds between atoms; indexed by particle index
			# structure: [atom,bonded-to,bond-length,gamma(bond-strength)]
		chrg			  = [] 
		positions 		  = []
		velocities 		  = []
		angles		 	 = {}	# dict tracking bond angles; indexed by particle index

		# iterate through and initialise values for all children (atoms in molecule)
		for i in range(self.n):
			# temporarily, this is set to extend the bond length array with duplicate values
			if i > len(bl)-1:
				bl.append(bl[i-1])
			if i > 0:
				# increment temporary position variable by bond length to give smoother load-in
				position = position + np.array([1.0/bl[i]**3,1.0/bl[i]**3,0]) 
				# in future, this will be replaced with / followed by Maxwell-Boltzmann to match initial velocities to Temp
			positions.append(position)
			#print(position)
			# initialise atom type, constants, absolute and relative positions, velocities, masses, charges, bonds, and angles
			childIndex 	= i + self.indx 
			child 		= self.children[i]
			#print(child)
			child.setIndx(childIndex)
			# the set of atoms in the arrays start at the molecule index; the molecule indexes the first atom in the array
			child.world_position = position
			child.velocity 		 = velocity
			velocities.append(velocity)
			self.mm.append(child.mass)
			chrg.append(child.charge)
			# assign the atom's parentMol variable to the new molecule
			child.parentMol = self
			# initialise bonds
			if i > 0:
				covbonds[childIndex] = [i-1,i,bl[i],Kr]
				if i < self.n:
					angles[child.tpindex] = [i-1,i,i+1,th0,Kth]

		# load lists into arrays
		RescalePosArray(positions)
		RescaleVelArray(velocities)
		RescaleMassArray(self.mm)
		self.charges 		= np.array(chrg)
		self.bl 			= np.array(bl)

		# bonds only required for molecules with more than 1 member
		if self.n > 1:
			self.covbonds 	= covbonds # alternatively, np.array(covbonds.values())
			# bond angles form between two bonds, so are not required if there are less than that
			if self.n > 2:
				self.angles = angles # alternatively, np.array(angles.values())
			else: # if there is only one bond in the molecule
				self.angles = [] # alternatively, np.array() 
		else: # if there are no bonds in the molecule
			self.covbonds = [] # alternatively, np.array() 

	def __array__():
		sumarray = np.array([self.children,self.mm,self.charges])

	def update(self):
		#global atomPositions
		#global atomVelocities
		#global atomAccel
		for i in range(self.n):
			self.children[i].world_position = atomPositions[self.indx+i]
			self.children[i].velocity 		= atomVelocities[self.indx+i]

	def positions(self):
		return [self.children[i].world_position for i in range(self.n)]

	def velocities(self):
		return [self.children[i].velocity for i in range(self.n)]


def RescalePosArray(*posits):
	global atomPositions
	# 	rescale atomPositions array
	if len(atomPositions) > 0:
		atmPos = atomPositions.tolist()
	else: atmPos = []
	# for each position array passed
	for pos in posits:
		print(pos)
		# iterate through the positions in the array passed and append them to the atomPositions array
		for po in pos:
			atmPos.append(po)
	atomPositions = np.array(atmPos)
	n = len(atomPositions)
	#print(atomPositions)

def RescaleVelArray(*velocs):
	global atomVelocities
	# 	rescale velocities array
	if len(atomVelocities) > 0:
		atmVel = atomVelocities.tolist()
	else: atmVel = []
	# for each velocity array passed
	for vel in velocs:
		# iterate through the velocities in the array passed and append them to the atomVelocities array
		for ve in vel:
			atmVel.append(ve)
	atomVelocities = np.array(atmVel)

def RescaleMassArray(*newmas):
	global masses
	# 	rescale velocities array
	if len(masses) > 0:
		masss = masses.tolist()
	else: masss = []
	for mas in newmas:
		for ma in mas:
			masss.append(ma)
	masses = np.array(masss)

def RescaleSigmArray(*newsig):
	global sigm
	# 	rescale velocities array
	if len(sigm) > 0:
		si = sigm.tolist()
	else: si = []
	for sii in newsig:
		for sigi in sii:
			si.append(sigi)
	sigm = np.array(si)

def RescaleEpsiArray(*neweps):
	global epsi
	# 	rescale velocities array
	if len(epsi) > 0:
		ep = epsi.tolist()
	else: ep = []
	for epp in neweps:
		for eppi in epp:
			ep.append(eppi)
	epsi = np.array(ep)

def CreateMolecule(location, *mols, player=0, bl=[1], velocity=[np.zeros(D)]):
	global molecules
	global n

	atfact = ats.AtomFactory() # factory method for generating atoms

	# update location based on player info
	if (player == 0 and location == None): location = Vec3(-30+randint(-5,5),-12+randint(-1,1),-20+randint(-5,5))
	elif (player == 1 and location == None): location = Vec3(30+randint(-5,5),-12+randint(-1,1),20+randint(-5,5))

	# run the creation procedure for every molecule passed in *mols tuple
	# LOOP PER MOL PASSED TO CONSTRUCT
	for mol in mols:
		nam = mol
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
				print("created " + newt.nam)
		# initialise arrays for sigma, epsilon, and type index values for the new atoms in the molecule
		newsig 	= []
		neweps 	= []
		tpis 	= []
		
		# IS THE SIMULATION EMPTY?
		if n > 0:
			# if the simulation contains molecules, generate new LJ parameters using the Lorentz-Berthelot combining rules
			for i in range(len(newatts)):
				# 0 generate type name
				name = mol + "_" + str(newatts[i].nam)
				# 1 initialise new index
				tpindex  = len(tp)
				# 2 search typearray for atom type
				for typ in range(len(tp)):
					if tp[typ] == name:
						tpindex = typ
				# type not found loaded: load new type, leave atom.unique = True
				if tpindex == len(tp):
					# 3 add type to array when loaded
					tp.append(name)
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
					for si in sigm: si.tolist().append((si[0]+newatts[i].sig)/2) # Lorentz: (sig_ii+sig_jj)/2 
					for ep in epsi: ep.tolist().append((ep[0]*newatts[i].eps)**.5) # Berthelot: sqrt(eps_ii*eps_jj) 
				else: # type found in tp array: save type index to atom object for easy lookup and log that atom is only one of its type
					newatts[i].tpindex = tpindex
					tpis.append(tpindex)
					# tell engine that atom is not unique in the simulation
					newatts[i].setUnique(False)
		else: # if this is the first molecule added to the sim, initialise arrays to molecule values and leave atom.unique = True
			for i in range(len(newatts)):
				# 0 generate type name
				name = mol + "_" + str(newatts[i].nam)
				# 1 initialise new index
				tpindex  = len(tp)
				# 2 search typearray for atom type
				for typ in range(len(tp)):
					if tp[typ] == name:
						tpindex = typ
				# type not found loaded: load new type, leave atom.unique = True
				if tpindex == len(tp):
					# 1 save indices
					newatts[i].tpindex = tpindex
					tp.append(name)
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
		mol = Molecule(nam, location, len(molecules), newatts, tpis, velocity)
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

# Gradient of Lennard-Jones potential on one body
def dLJp(r,mol,atid,i,sigl,epsl,bdln):
	#atm = mol.children[atid]
	pvi = mol.children[atid].indx # position and vector array index of atom

	print(mol.children[atid].nam)
	print(mol.children[atid].world_position)
	print(mol.children[atid].velocity)
	print(mol.children[atid].sig)
	print(mol.children[atid].tpindex)
	sg = np.array([sigl[i]])
	# if this is the only atom of its type, remove its sigma value to create the temporary array
	if mol.children[atid].isUnique():
		#print(sigl)
		#print("sigl[i]")
		#print(sigl[i])
		sg=np.delete(sigl,0)
	
	ep = np.array([epsl[i]])
	# remove epsilon value for self if it is unique
	if mol.children[atid].isUnique(): ep=np.delete(ep,0) 
	
	# calculate distance vector then absolute distance in each dimension
	drv=r-r[pvi] # distance vector for every point
	#print ("DRV:")print(drv)
	for j in range(mol.n): # remove position vectors of all atoms in the molecule from drv
		drv= np.delete(drv,mol.indx,0) # n.b. tge index of the removed value doesn't change between iterations due to the resizing of the array
	#print ("DRV after molecule removal:")print(drv)
	dr=[np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) for a in drv] #absolute distance of that lad

	# truncate to only calculate only for nearby (less than 4 fm away) atoms. May be replaced by partitioning
	if ((dr[i] > 4) for i in dr):
		dLJP=0.0
	else: # calculate dLJP (8-14 from 6-12)
		#print("ep: " + str(ep))
		#print("dr: " + str(dr))
		# TODO
		# 	figure out how the hell to look up the correct sigma and epsilon values for the target interactions
		r8 = ep*(sg**6)*(1.0/np.array(dr))**8 # IF THERE'S A DIVIDE BY ZERO, you're calculating particle self-interactions
		r14=2.0*ep*(sg**12)*(1.0/np.array(dr))**14
		r8v =np.transpose(np.transpose(drv)*r8)
		r14v=np.transpose(np.transpose(drv)*r14)
		r8vs =np.sum(r8v,axis=0)
		r14vs=np.sum(r14v,axis=0)
		dLJP=24.0*(r14vs-r8vs)

	return dLJP

#bond length potential
def BEpot(r,bnds):
	bps=np.zeros(n)
	print(bnds)
	for i in range(n): #loop over all particles
		for j in range(len(bnds)): #check all bonds to see if particle i is bonded
			if bnds[j][0]==i or bnds[j][1]==i:
				print("")
				if bnds[j][0]==i: #find particle bonded to i
					ii=int(bnds[j][1])
				else:
					ii=int(bnds[j][0])
				dr0=bnds[j][2]
				e0 =bnds[j][3]
				dr=r[i]-r[ii]
				dr2=dr*dr
				adr2=sum(dr2)
				adr=np.sqrt(adr2)
				BE=e0*(adr-dr0)**2
				bps[i]+=BE
	return bps

#gradient of bond length potential (negative force)
def dBEpot(r,bnds):
	bps=np.zeros([n,3])
	for i in range(n): #loop over all particles
		for j in range(len(bnds)): #check all bonds to see if particle i is bonded
			if bnds[j][0]==i or bnds[j][1]==i:
				if bnds[j][0]==i: #find particle bonded to i
					ii=int(bnds[j][1])
				else:
					ii=int(bnds[j][0])
				dr0=bnds[j][2]
				e0 =bnds[j][3]
				dr=r[i]-r[ii]
				dr2=dr*dr
				adr2=sum(dr2)
				adr=np.sqrt(adr2)
				dBE=2.0*e0*(adr-dr0)*dr/adr
				bps[i]+=dBE
	return bps

#gradient of bond angle potential (negative force)
def dBA(r,angs):
	aps=np.zeros([n,3])
	for i in range(n): #loop over all particles
		for j in range(len(angs)): #check all bonds to see if particle i is bonded
			a1=int(angs[j][0])
			a2=int(angs[j][1])
			a3=int(angs[j][2])
			if i==a1 or i==a2 or i==a3:
				th00=angs[j][3] #equilibrium angle
				e0 =angs[j][4] #bending modulus
				if i==a1 or i==a2:
					r1=r[a1]-r[a2] #bond vector 1 (form middle atom to atom 1)
					r2=r[a3]-r[a2] #bond vector 2 (middle atom to atom 2)
				else:
					r1=r[a3]-r[a2] #bond vector 1 (form middle atom to atom 1)
					r2=r[a1]-r[a2] #bond vector 2 (middle atom to atom 2)
				ar1=np.sqrt(sum(r1*r1)) #lengths of bonds
				ar2=np.sqrt(sum(r2*r2))
				dot=sum(r1*r2) #r1 dot r2
				ndot=dot/(ar1*ar2) #normalize dot product by vector lengths i.e. get the cos of angle
				th=math.acos(ndot) #bond angle, theta
				dUdth=-2.0*e0*(th-th00) #-dU/dtheta
				if a1==i or a3==i:
					numerator=(r2/(ar1*ar2))-(dot/(ar1*ar1*ar1*ar2*2.0))
					denominator=np.sqrt(1.0-ndot*ndot)
					dUdr=dUdth*numerator/denominator
					aps[i]+=dUdr
				if i==a2:
					denominator=np.sqrt(1.0-ndot*ndot)
					n1=-(r2+r1)
					n2=dot*r1/(ar1*ar1)
					n3=dot*r2/(ar2*ar2)
					numerator=(n1+n2+n3)/(ar1*ar2)
					dUdr=dUdth*numerator/denominator#
					aps[i]+=dUdr
	return aps

#derivative of coulomb potential (negative force)
def coul(r,i,q0,mols):
	qs=1.0*np.array([ml.chrg for ml in molecules])
	for j in range(n):
		if qs[i]==mols[j]:
			qs[j]=0.0
	qs=np.delete(qs,i)
	drv=r-r[i] #distance in each dimension
	drv=np.delete(drv,i,0) #remove ith element (no self LJ interactions)
	dr=[np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) for a in drv] #absolute distance
	r3=q0*qs*kc*((1.0/np.array(dr))**3.0) # Coulomb's law
	FF =np.transpose(np.transpose(drv)*r3) # transpose the distance array, multiply by force, transpose back
	Fs=np.sum(FF,axis=0) # sum the value of axis 0 of the transposed & multiplied array
	print(Fs)
	return Fs
	#global D
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
	j = 0
	a = []
	if n > 0:
		for mol in molecules:
			# TODO calculate values in parallel
			lj = -np.array([dLJp(atomPositions,mol,atm,mol.children[atm].tpindex,sigm,epsi,mol.bl[atm]) for atm in range(mol.n)])
			bep= -dBEpot(atomPositions,mol.covbonds) #Bonds
			ba = -dBA(atomPositions,mol.covbonds) #Bond angles
			ch = -np.array([coul(atomPositions,mol.children[atm].indx,mol.chrgs[atm],molecules) for atm in range(mol.n)]) #Coulomb
			F= ((lj + bep) + ba) + ch
			a[j]=np.transpose(np.transpose(F)/mol.mm) #Force->acceleration
			j += 1
		a = np.array(a)

		# update velocity !!!					!!! POSSIBLE BUG
		atomVelocities = atomVelocities + time.dt * a # I'm not sure if using this new 2D accelleration array will work in this equation as-is
		atomVelocities = rescaleT(atomVelocities,Temp0) # scale to temperature
		atomPositions = atomPositions + atomVelocities * time.dt  # should be average of atomVelocities over time.dt
		print("updated " + mol.nam)
		# split up calculation if dt gets too large to prevent accumulating errors

		# boundaries
		if BC == True:
			atomPositions = atomPositions%L
			atomVelocities = 1.0 * atomVelocities
		else:		
			atomPositions, atomVelocities = reflectBC(atomPositions,atomVelocities)
		#total_T_display.text = calculateTemperature(atomVelocities,n)
		#collisions_display.text = str(collisions_count)
		#mouseloc_display.text = str(mouse.collisions)
