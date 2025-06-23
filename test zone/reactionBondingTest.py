from ursina import *
import numpy as np
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
DirectionalLight(parent=pivot, y=1.5, z=3, shadows=True, rotation=(65, -15, 45))

# some physical constants
kb  = 0.8314459920816467 # Boltzmann
NA  = 6.0221409e+26 #Avogardos constant x 1000 (g->kg)
ech = 1.60217662E-19 #electron charge in coulombs

# constant deltatime
setdt = 0.005

# game entities with initial position vectors
firstAtom = Entity(model='sphere', scale=0.1, position=np.array([-2,0,0]), color=color.red)
secondAtom = Entity(model='sphere', scale=0.1, position=np.array([2,0,0]), color=color.blue)

# derivatives of motion for calculations
veli = np.array([0.0,0.0,0.0])
velj = np.array([0.0,0.0,0.0])
acci = np.array([0.0,0.0,0.0])
accj = np.array([0.0,0.0,0.0])

# physical values
chrgi = 0.41
chrgj = 0.41
massi = 1.0
massj = 1.0

# dLJP for Hydrogen
def dLJP(parti,partj):
	ep = 3.24
	sg = 0.98

	drv  = partj.position - parti.position
	dr   = np.sqrt(drv[0]**2+drv[1]**2+drv[2]**2)
	r8vs = np.sum(np.transpose(np.transpose(drv) * ep*(sg**6) * (1.0/np.array(dr))**8),axis=0)
	r14vs= np.sum(np.transpose(np.transpose(drv) * 2.0*ep*(sg**12) * (1.0/np.array(dr))**14),axis=0)
	return 24.0*(r14vs-r8vs)

def coul(parti,partj,qi,qj):
	kc = 8.9875517923E9*NA*1E30*ech*ech/1E24 #electrostatic constant in Daltons, electron charges, picosecond, angstrom units

	drv = partj.position - parti.position
	dr  = np.sqrt(drv[0]**2 + drv[1]**2 + drv[2]**2) 		 # absolute distance of those lads
	r3  = kc * qi * qj * ((1.0/dr)**3.0) 					 # Coulomb's law
	return np.sum(np.transpose(np.transpose(drv)*r3),axis=0) # transpose the distance vector, multiply by force, transpose back

def update():
	global veli, velj, acci, accj

	# force calculations
	forcei = -np.array([dLJP(firstAtom,secondAtom) + coul(firstAtom,secondAtom,chrgi,chrgj)])
	forcej = -np.array([dLJP(secondAtom,firstAtom) + coul(secondAtom,firstAtom,chrgj,chrgi)])

	# a = F/m
	acci = np.transpose(np.transpose(forcei)/massi)
	accj = np.transpose(np.transpose(forcej)/massj)

	print(forcei)
	print(acci)
	print(veli)

	# euler integration of position derivatives to position
	veli = veli + (setdt * acci)
	velj = velj + (setdt * accj)
	firstAtom.position = firstAtom.position + (veli * setdt)
	secondAtom.position = firstAtom.position + (velj * setdt)

	camera.lookAt(firstAtom)

app.run()
