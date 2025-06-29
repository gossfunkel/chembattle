# basic physics sim with discretely defined entities (no array)
# this file is for testing different models of chemical bonding simulation
# 	primarily focussed on relationships between atoms/nuclei

#TODO :
# TEST 1: PEP
# set up a carbon nucleus with no electrons. Add electrons slowly so that they fill the potential well of the 
# nucleus quantum by quantum (filling orbital shells)
# TEST 2: HYBRIDISATION
# Bring two Hydrogen ions near a Sulfur ion such that the S2 orbital of the Sulfur will have to hybridise with 
# two p-orbitals to allow a bonding configuration. This raises the potential energy level of the two new sp orbitals.
# Then try Nitrogen and 3x Hydrogen for ammonia, and 
# TEST 3: PI BONDS
# Try the following with CO2 as well to see how it works on linear geometry;
# Construct ethylene (2H sigma bonded to an sp3 hybridised C bonded to another unit of the same, mirrored)
# Test out systems for electronic repulsion such that the p-orbital electrons localise to the pi-bond space
# nuclei attract electrons, electrons repel each other, pep defines the limits
# slow dt down so far; something like 6.0e-17, such that nuclei barely move and electrons wander slowly
# TEST 4: RESONANCE STRUCTURES
# Construct a benzene ring. Ensure that the electrons localise appropriately to the two rings. Note that there
# should now be a small EM field through the hole in the toroidal body of the molecule. 
# TEST 5: REDOX
# Can you get a Sodium to give away an electron to Chlorine, on the terms that they become ionically associated
# as atoms due to the new net charge on each atom? 

from ursina import *
import numpy as np
#from panda3d.core import Shader
#from chempy import Substance

# ursina boilerplate to set up window, camera, and lighting
WINDOW_WIDTH, WINDOW_HEIGHT = 1080, 600
app = Ursina(size=(WINDOW_WIDTH,WINDOW_HEIGHT))
window.vsync = False
window.title = "BONDING POTENTIAL TEST"
window.borderless = False
window.fullscreen = False
window.fps_counter.enabled = True
window.exit_button.enabled = False
window.color = color.black
EditorCamera()
sky = Sky()
pivot = Entity()
#DirectionalLight(parent=pivot, y=1.5, z=3, shadows=True, rotation=(65, -15, 45))
selectView = True

testShader = Shader(Shader.GLSL, vertex="testVertexShader.vert", fragment="testFragShader.frag")

# some physical constants
kb  = 0.8314459920816467 # Boltzmann
NA  = 6.0221409e+26 #Avogardos constant x 1000 (g->kg)
ech = 1.60217662E-19 #electron charge in coulombs
gravConst = 6.6743e-11 # G in m3/kg/s

# constant deltatime
setdt = 0.00001

# game entities with initial position vectors
firstAtom = Entity(model='sphere', scale=1., world_position=np.array([1,-10,35]), color=color.red, shader=testShader)
secondAtom = Entity(model='sphere', scale=1., world_position=np.array([14,-10,15]), color=color.blue, shader=testShader)
camera.world_position = Vec3(15,0,22)
camera.shader = testShader
#camera.rotation_x = 0

#visualPlane = Entity(model='plane',collider=None, world_position=np.array([10,-10,20]))

# derivatives of motion for calculations
veli = np.array([-50.0,0.0,-50.])
velj = np.array([50.0,0.0,50.])
#acci = np.array([0.0,0.0,0.0])
#accj = np.array([0.0,0.0,0.0])

# physical values
chrgi = 0.41
chrgj = 0.41
massi = 0.001
massj = 0.001

# dLJP for Hydrogen
def dLJP(parti,partj):
	ep = 3.24
	sg = 0.98
	#print(parti)
	drv  = partj - parti
	dr   = np.sqrt(drv[0]**2 + drv[1]**2 + drv[2]**2)
	#r8vs = np.sum(np.transpose(np.transpose(drv) * ep*(sg**6) * (1.0/np.array(dr))**8),axis=0)
	r8vs = np.transpose(np.transpose(drv*(ep*(sg**6) * (1.0/np.array(dr))**8)))
	r14vs= np.transpose(np.transpose(drv) * 2.0*ep*(sg**12) * (1.0/np.array(dr))**14)
	#r14vs = np.sum(drv*(2.0*ep*(sg**12) * (1.0/np.array(dr))**14))
	return np.array((r14vs-r8vs)*24.0)

def coul(parti,partj,qi=0.41,qj=0.41):
	kc = 8.9875517923E9*NA*1E30*ech*ech/1E24 #electrostatic constant in Daltons, electron charges, picosecond, angstrom units
	drv = partj - parti
	dr  = np.sqrt(drv[0]**2 + drv[1]**2 + drv[2]**2) 		 # absolute distance of those lads
	r3  = kc * qi * qj * (1.0/dr**3.0) 					 	 # Coulomb's law
	return np.transpose(np.transpose(drv)*r3)				 # transpose the distance vector, multiply by force, transpose back
	#return np.sum(drv*r3)

def grav(parti,partj):
	drv = partj - parti
	dr  = np.sqrt(drv[0]**2 + drv[1]**2 + drv[2]**2) 		 # absolute distance of those lads
	return gravConst * massi * massj / np.transpose(np.transpose(dr)**2)

def ode(rx, ry, v0, dt):
	#forces = np.zeros(2,3)
	force = -dLJP(rx, ry) + coul(rx, ry) - grav(rx,ry)
	#print(force)
	#forcej = -np.array([dLJp(rj, ri) - coul(rj, ri)])
	#a = np.transpose(np.transpose(forces)/masses]) # find r''(dt) = a(r,q,m)
	a = np.transpose(np.transpose(force) / massi) 	# !!!!!!!!!!!!!!!!! hacky - only works while masses are equal
	v = v0 + a * dt 								# find r'(dt) = v(dt) from a 
	rx = rx + v * dt 								# find r(dt)
	ry + ry - v * dt
	return a,v,rx,ry 									# return dv/dt and v for rk4, return r for newton
	
def rk4(posi, posj, velo, dt):
    k1,velo,posi, posj = ode(posi, posj, velo, 0)
    k2,velo,posi, posj = ode(posi, posj, velo + k1*dt/2, dt/2)
    k3,velo,posi, posj = ode(posi, posj, velo + k2*dt/2, dt/2)
    k4,velo,posi, posj = ode(posi, posj, velo + k3*dt, dt)
    return posi, velo + dt/6*(k1 + 2*k2 + 2*k3 + k4)

def input(self):
	global selectView
	if held_keys['space']:
		selectView = not selectView

def update():
	global veli, velj
	firstAtom.position, veli = rk4(firstAtom.position, secondAtom.position, veli, setdt)
	secondAtom.position, velj = rk4(secondAtom.position, firstAtom.position, velj, setdt)
	print(veli)
	# force calculations
	#forcei = -np.array([dLJP(firstAtom,secondAtom) + coul(firstAtom,secondAtom,chrgi,chrgj)])
	#forcej = -np.array([dLJP(secondAtom,firstAtom) + coul(secondAtom,firstAtom,chrgj,chrgi)])
	# a = F/m
	#acci = np.transpose(np.transpose(forcei)/massi)
	#accj = np.transpose(np.transpose(forcej)/massj)
	# euler integration of position derivatives to position
	#veli = veli + (setdt * acci)
	#velj = velj + (setdt * accj)
	#firstAtom.position = firstAtom.position + (veli * setdt)
	#secondAtom.position = firstAtom.position + (velj * setdt)
	
	if selectView: camera.lookAt(firstAtom)
	else: camera.lookAt(secondAtom)

app.run()
