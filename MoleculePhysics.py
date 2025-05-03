from ursina import *
#from random import uniform
from random import uniform
import numpy as np

#molecules = np.empty([2,0], np.dtype(object))
molecules = []
dragstate = False

def AddMolecule(Molecule, player=0):
	# add to correct axis of array for which player is adding the molecule
	# default to player 1 (local) 
	molecules.append(Molecule)
	print(Molecule.children)

def GetMolecules(player=0):
	return molecules
	# return a 1D array of all the molecules for a player
	# default to player 1 (local)
	#if (molecules.shape[player] > 0):
	#	return molecules[player, :]
	#else: return None

def CreateMolecule(location, velocity, atom, player=0):
	# initiate a molecule with parameters passed
	newIon = Molecule(location, velocity, atom)
	# assign the atom's parentMol variable to the new molecule
	atom.parentMol = newIon
	# add the molecule to the numpy array of all molecules
	AddMolecule(newIon, player)
	# return the molecule object for any further assignments etc
	return newIon

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

def React(mola, molb):
	mola.velocity = -mola.velocity
	molb.velocity = -molb.velocity

class Molecule(Draggable):
	#def drag(self):
	#	self.dragging = True
	#	for cld in self.children:
	#		cld.dragging = True

	def __init__(self, position, velocity, *members):
		#if (len(members) == 0):
		#	return None
		super().__init__(position=position)
		#self.npposition = np.array([position])
		#self.world_position  = Vec3(position[0],position[1],position[2])
		#print(position)
		self.velocity     = velocity
		#self.dragging     = False
		self.children     = members
		self.visible_self = False
		#self.collider    = 'sphere'
		self.mass 	      = 0
		self.temp 		  = 0
		self.rad 		  = 0
		for mem in members:
			# mass of molecule is the sum of the mass of its' constituent atoms' masses
			self.velocity += mem.velocity
			self.mass    += mem.mass
			# temporary radius of molecule is the sum of its' constituent atoms' radii
			self.rad     += mem.rad
			# temperature of molecule is the average of its' constituent atoms' temperatures
			self.temp    += mem.temp 
		self.temp         = self.temp/len(members)
		self.wobble       = (self.temp*self.rad)
		#self.emforce  = np.zeros(3)
		#self.prforce  = np.zeros(3)

	def update(self):
		self.movement()

	def movement(self):
		# Coulomb's law: f = k*q1*q2/r**
		#k = 0.000000008988 # Nm**/C**
		# however, the Lennard-Jones potential is better worked out via:
		# Van der Waals' attractive force: F a r**6
		# Pauli's repulsive force: F a 1/x**12
		#for ch in charges:
		#	mol.emforce += k * mol.charge * ch.charge / pygame.math.distance(mol.position, ch.position)**
		#mol.accel = mol.mass * (mol.emforce + mol.prforce)
		#mol.velocity *= mol.accel
		#print(mol.position)
		#print(mol.velocity[0])
		for chi in self.children:
			self.velocity += chi.velocity
		#self.velocity += chi in self.children

		self.world_position += self.velocity * time.dt
		#rotationfactor = Vec3((np.pi * time.dt),(np.pi * time.dt),(np.pi * time.dt))
		#self.rotate(rotationfactor)
		self.wobblevec = Vec3(uniform(-self.wobble,self.wobble),uniform(-self.wobble,self.wobble),uniform(-self.wobble,self.wobble))
		#self.position = SlideTo(self.position + self.velocity, self.position, 15)
		self.position = SlideTo(self.position + self.wobblevec, self.position, 15)
		#print(self.position)
		#print(self.velocity)
		#print(self.wobble)

#	def react(self, mol):
		# pseudocode:
		#if (self.satisfiesOctetRule && mol.satisfiesOctetRule):
		#	repel at close range
		#else:
		# 	if ((self.charge + mol.charge) == 0):
		#		CalculateMolecularReaction(self, mol)
		#	else:
		#		if (self.charge > mol.charge):
		#			electrophile = self
		#			nucleophile = mol
		#		else:
		#			electrophile = self
		#			nucleophile = mol
		#		CalculateIonicMolecularReaction(electrophile, nucleophile)

# ''' from ursina api ref for CubicBezier
# Curves from Ursina are used by Entity when animating, like this:

# e = Entity()
# e.animate_y(1, curve=curve.in_expo)

# e2 = Entity(x=1.5)
# e2.animate_y(1, curve=curve.CubicBezier(0,.7,1,.3))
# '''