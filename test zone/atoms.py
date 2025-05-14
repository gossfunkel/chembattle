from ursina import *
import chempy as chem
import numpy as np
import AdvancedMoleculePhysics as amp
from os.path import join
from random import randint

#shells = ['s1','s2','p1','p2','p3','s3','s4','d1','d2','d3','d4','p4','p5','p6']
#textures = {
#	'Hydrogen' : 
#}

def CreateAtom(atom, player=0):
	ion   = 0
	veloc = Vec3(0.6,0,0.6)
	if (player == 0):
		player = 1
		locat = Vec3(-30,-12,-20)
	elif (player == 1):
		player = -1
		locat = Vec3(30,-12,20)

	if (atom == 'Hydrogen'):
		newat = Hydrogen(locat, ion, veloc)
	elif (atom == 'Carbon'):
		locat += Vec3(0,0,player*2)
		newat = Carbon(locat, ion, veloc)
		#print(locat)
	elif (atom == 'Nitrogen'):
		locat += Vec3(0,0,player*4)
		newat = Nitrogen(locat, ion, veloc)
	elif (atom == 'Oxygen'):
		locat += Vec3(0,0,player*6)
		newat = Oxygen(locat, ion, veloc)
	else: return None
	
	newat.parent = amp.CreateMolecule(newat.name, locat, newat.sig, newat.eps, newat.bl, newat, player)
	print("adding " + atom)
	return newat

class Electron(Entity):
 	def __init__(self, position, binding, eshell=np.array([1.0,0.0,0.0]), spin=0):
 		super().__init__(collides=False,model='sphere',scale=(0.1,0.1,0.1))
 		# spin opposite electrons are on the opposite sides for now
 		if spin == 0: self.eshell = eshell
	 	else: self.eshell = -eshell
 		#join('res',shells[eshell] + '.png')  # pick correct image using dictionary of energy shells
 		self.world_position = position + eshell
 		self.binding = binding
 		self.t = -np.pi
 		self.spin = 0.5 - spin
 		#self.laz = self.eshell[0]
 		#self.n = self.eshell[1]
 		#self.mag = self.ehell[3]

 	def update(self):
 		# temporarily, the electrons shall orbit the atoms
 		# omicron = 2pi/T where omicron is the angular velocity and T is the period
 		# angle theta = 2pi*t/T where t is the time elapsed
 		# acceleration due to change in direction = v**/r = omicron**r
 		# so we can calculate the new position as the change of angle of rotation in the electron's position vector's directional component
 		# then we can use a rotation matrix transform i found on stackexchange to change the position modification vector eshell:
 		# |  x*cos(theta) + z*sin(theta) |
 		# |              y				 |
 		# | -x*sin(theta) + z*cos(theta) |
  		# variables for orbit
 		self.t += 0.2 # period of rotation
 		theta = 2*np.pi*time.dt/2 # period of 2
 		# turn theta into a new displacement vector
 		self.eshell = Vec3(self.eshell[0]*np.cos(theta) + self.eshell[2]*np.sin(theta),self.eshell[1],-self.eshell[0]*np.sin(theta) + self.eshell[2]*np.cos(theta))
 		
 		# temporary 'total potential energy' of electron "calculated" from quantum numbers
 		#nrg = 1 * self.eshell[0] + 1 * self.eshell[1] + 2 * self.eshell[2]
 		# more temporary integer used as energy value
 		#nrg = self.eshell[0]
		
 		self.world_position = self.parent.world_position + self.eshell
 		# note: spin opposition in space is coded for in init: if statement reverses eshell direction depending on spin value
			#self.x = self.parent.x + theta + 180
 			#self.z = self.parent.z + theta 
			#np.cos(self.t)*nrg*time.dt
 			#np.sin(self.t)*nrg*time.dt
 			#self.x = np.cos(self.parent.x)*nrg*time.dt
 			#self.z = np.sin(self.parent.z)*nrg*time.dt

class Atom(Entity):
	def __init__(self, position, scale=(1,1,1), ionisation=0, velocity=np.zeros(3), uri="default.png", temp=0.0, electrons=[]):
		#super().__init__(billboard=True)
		super().__init__(model='sphere',visible=True)
		self.world_position   = position
		self.ionisation = ionisation
		self.velocity   = velocity
		self.rad        = self.scale[0]
		self.temp       = temp
		#self.model 		= 'atoms' blender is being stroppy so not yet
		self.texture    = join('textures',uri)
 		self.sig 		= 1 # TODO Van der Waals radius
 		self.eps 		= 1 # TODO Energy well depth
 		self.bl 		= 1 # TODO Spring potential equilibrium radius
		
		# may have to switch this out and calculate based on exact charge of proton in eV
		self.charge = self.atomNum + self.ionisation
		nextSpin = False
		for ele in range(self.charge):
			if not nextSpin:
				nextSpin = True
				e = Electron(position,self,Vec3(ele+1,0,0),0)
			else:
				nextSpin = False
				e = Electron(position,self,Vec3(ele,0,0),1)
			self.children.append(e)
			e.parent = self
		#self.mass = mass
		#self.hybridisation = hybridisation

class Hydrogen(Atom):
	def __init__(self, position, ionisation=0, velocity=np.zeros(3), temp=0.04):
		self.atomNum = 1
		self.uri     = 'atomHTrans.png'
		self.mass    = 1.008
		self.size    = (2.5,2.5,2.5)
		super().__init__(position, self.size, ionisation, velocity, self.uri, temp)
		self.color  = '#028cad'
		self.sig	= 1.2

class Lithium(Atom):
	def __init__(self, position, ionisation=0, velocity=np.zeros(3), temp=0.04):
		self.atomNum = 3
		self.uri     = 'atomLiTrans.png'
		self.mass    = 6.946
		self.size    = (2.8,2.8,2.8)
		super().__init__(position, self.size, ionisation, velocity, self.uri, temp)
		self.color  = '#ffbbb0'
		self.sig	= []

class Carbon(Atom):
	def __init__(self, position, ionisation=0, velocity=np.zeros(3), temp=0.04):
		self.atomNum = 6
		self.uri     = 'atomCTrans.png'
		self.mass    = 12.011
		self.size    = (3.8,3.8,3.8)
		super().__init__(position, self.size, ionisation, velocity, self.uri, temp)
		self.color  = '#c8f900'

class Nitrogen(Atom):
	def __init__(self, position, ionisation=0, velocity=np.zeros(3), temp=0.04):
		self.atomNum = 7
		self.uri     = 'atomNTrans.png'
		self.mass    = 14.007
		self.size    = (4.0,4.0,4.0)
		super().__init__(position, self.size, ionisation, velocity, self.uri, temp)
		self.color  = '#10ff06'

class Oxygen(Atom):
	def __init__(self, position, ionisation=0, velocity=np.zeros(3), temp=0.04):
		self.atomNum = 8
		self.uri     = 'atomOTrans.png'
		self.mass    = 15.999
		self.size    = (4.2,4.2,4.2)
		super().__init__(position, self.size, ionisation, velocity, self.uri, temp)
		self.color  = '#ff1005'

class Phosphorus(Atom):
	def __init__(self, position, ionisation=0, velocity=np.zeros(3), temp=0.04):
		self.atomNum = 15
		self.uri     = 'atomPTrans.png'
		self.mass    = 30.973
		self.size    = (4.8,4.8,4.8)
		super().__init__(position, self.size, ionisation, velocity, self.uri, temp)
		self.color  = '#ff0000'

class Sulfur(Atom):
	def __init__(self, position, ionisation=0, velocity=np.zeros(3), temp=0.04):
		self.atomNum = 16
		self.uri     = 'atomSTrans.png'
		self.mass    = 32.062
		self.size    = (5.2,5.2,5.2)
		super().__init__(position, self.size, ionisation, velocity, self.uri, temp)
		self.color  = '#00ff00'

class Chlorine(Atom):
	def __init__(self, position, ionisation=0, velocity=np.zeros(3), temp=0.04):
		self.atomNum = 17
		self.uri     = 'atomClTrans.png'
		self.mass    = 35.451
		self.size    = (5.5,5.5,5.5)
		super().__init__(position, self.size, ionisation, velocity, self.uri, temp)
		self.color  = '#ffff00'

#class Silicon(Atom):
#	def __init__(self, position, )