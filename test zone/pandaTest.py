from direct.showbase.ShowBase import ShowBase
from panda3d.core import Shader, Geom, GeomNode, GeomTriangles, GeomVertexWriter, loadPrcFileData
from panda3d.core import GeomVertexFormat, GeomVertexData, TransparencyAttrib, NodePath
from math import sin, cos
import numpy as np

# def createColouredRect(x, z, width, height, colors=None):
	# 	format = GeomVertexFormat.getV3c4()
	# 	vdata   = GeomVertexData('square', format, Geom.UHDynamic)

	# 	vertex  = GeomVertexWriter(vdata, 'vertex')
	# 	color   = GeomVertexWriter(vdata, 'color')

	# 	vertex.addData3(x, 0, z)
	# 	vertex.addData3(x+width, 0, z)
	# 	vertex.addData3(x+width, 0, z+height)
	# 	vertex.addData3(x, 0, z+height)

	# 	if colors:
	# 		if len(colors) < 4:
	# 			colors = (colors[0],colors[0],colors[0],1.0)
	# 		color.addData4f(colors)
	# 		color.addData4f(colors)
	# 		color.addData4f(colors)
	# 		color.addData4f(colors)
	# 	else:
	# 		color.addData4f(1.0,0.0,0.0,0.0)
	# 		color.addData4f(0.0,1.0,0.0,0.0)
	# 		color.addData4f(0.0,0.0,1.0,0.0)
	# 		color.addData4f(0.0,0.0,0.0,1.0)

	# 	tris = GeomTriangles(Geom.UHDynamic)
	# 	tris.addVertices(0,1,2)
	# 	tris.addVertices(2,3,0)

	# 	square = Geom(vdata)
	# 	square.addPrimitive(tris)
	# 	return square
  
loadPrcFileData('', 'win-size 1000 600') 
loadPrcFileData('', 'show-frame-rate-meter 1')
loadPrcFileData('', 'hardware-animated-vertices true')
loadPrcFileData('', 'basic-shaders-only false')

	## how many spheres
sphereNum = 25 

	### physical values
chrg = np.array([0.41 for _ in range(sphereNum)])
mass = np.array([0.001 for _ in range(sphereNum)])
	### some physical constants
kb  = 0.8314459920816467 # Boltzmann
NA  = 6.0221409e+26 #Avogardos constant x 1000 (g->kg)
ech = 1.60217662E-19 #electron charge in coulombs
kc = 8.9875517923E9*NA*1E30*ech*ech/1E24 #electrostatic constant in Daltons, electron charges, picosecond, angstrom units
#gravConst = 6.6743e-11 # G in m3/kg/s
	## constant deltatime - n.b. for scientific accuracy, we would want this to have a resolution of 
setdt = 0.0000003
	## derivatives of motion for calculations
vel = np.array([[i,i,0.0] for i in range(sphereNum)])
pos = np.array([[(sin(i)*30)+30,-cos(i)*30.0,5.0+i/2.0] for i in range(sphereNum)])
#ep = np.array([3.24 for _ in range(sphereNum)])
#sg = np.array([0.98 for _ in range(sphereNum)])

# dLJP for Hydrogen
def dLJP(r,i):
	ep = 3.24
	sg = 0.98
	#print(parti)
	drv  = r - r[i]
	drv=np.delete(drv,i,0) #remove ith element (no self LJ interactions)
	dr=[np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) for a in drv] #absolute distance of that lad
	#r8vs = np.sum(np.transpose(np.transpose(drv) * ep*(sg**6) * (1.0/np.array(dr))**8),axis=0)
	dLJP=0.0
	if ((dr[i] < 4) for i in dr): # calculate dLJP (8-14 from 6-12)
		r8vs = np.sum(np.transpose(np.transpose(drv)*(ep*(sg**6) * (1.0/np.array(dr))**8)),axis=0)
		r14vs= np.sum(np.transpose(np.transpose(drv) * 2.0*ep*(sg**12) * (1.0/np.array(dr))**14),axis=0)
		dLJP = (r14vs-r8vs)*24.0
	return dLJP

def coul(r,i,qi=0.41,qj=0.41):
	drv = r - r[i]
	drv=np.delete(drv,i,0) #remove ith element (no self EM interactions)
	dr=[np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) for a in drv] #absolute distance of that lad
	r3  = kc * qi * qj * ((1.0/np.array(dr))**3.0) 					 	 	 # Coulomb's law
	return np.sum(np.transpose(np.transpose(drv)*r3), axis=0)				 	 # transpose the distance vector, multiply by force, transpose back

def ode(r, v, dt):
	forces = np.zeros((sphereNum,3))
	#print(force)
	#forcej = -np.array([dLJp(rj, ri) - coul(rj, ri)])
	#a = np.transpose(np.transpose(forces)/masses]) # find r''(dt) = a(r,q,m)
	for i in range(sphereNum):
		forces[i] = -dLJP(r,i) + coul(r,i) #- grav(rx,ry)
		a = np.transpose(np.transpose(forces[i]) / mass[i]) 	# !!!!!!!!!!!!!!!!! hacky - only works while masses are equal
		v[i] = v[i] + a * dt 								# find r'(dt) = v(dt) from a 
		r[i] = r[i] + v[i] * dt 							# find r(dt)
	return a,v,r										# return dv/dt and v for rk4, return r for newton
	
def rk4(pos, vel, dt):
    k1,vel,pos = ode(pos, vel, 0)
    k2,vel,pos = ode(pos, vel + k1*dt/2, dt/2)
    k3,vel,pos = ode(pos, vel + k2*dt/2, dt/2)
    k4,vel,pos = ode(pos, vel + k3*dt, dt)
    return pos, vel + dt/6*(k1 + 2*k2 + 2*k3 + k4)

class TestBase(ShowBase):
	def __init__(self):
		ShowBase.__init__(self)
		# get nodepath for environment model, then reparent to render node
		#self.scene = self.loader.loadModel("models/environment")
		#self.scene.reparentTo(self.render)
		# Apply scale and position transforms on the model.
		#self.scene.setScale(0.25)
		#self.scene.setPos(0, 42, 0)

		self.cam.setPos(-5, 0, 30)

		#newSquare = createColouredRect(0, 0, 100, 100)

		#geomNode = GeomNode('square')
		#geomNode.addGeom(newSquare)
		#self.render.attachNewNode(geomNode)

		# for some reason, the fps seems to cap at 33 whether I generate 5 spheres or a few hundred
			# should be hardware instanced like this - not 100% sure how panda3d manages this
			# i feel like i should also probably use panda3d's Intervals to update the motion in parallel;
					# see https://docs.panda3d.org/1.10/python/programming/intervals/index#intervals
					# and https://docs.panda3d.org/1.10/python/introduction/tutorial/using-intervals-to-move-the-panda
		self.spheres = np.zeros(sphereNum)
		self.thetas = np.zeros(sphereNum)
		self.sphere = self.loader.loadModel("../sphere")
		self.sphereNodep = NodePath('spheres') # placeholders can be attached to another node to spawn groups of instances
		for i in range(sphereNum):
			placeholder = self.sphereNodep.attachNewNode("Sphere-Placeholder")
			placeholder.setScale(0.5)
			placeholder.setPos(pos[i,0],pos[i,1],pos[i,2])
			self.sphere.instanceTo(placeholder)
			self.thetas[i] = 0.0
		placeholder = render.attachNewNode("SphereGroup")
		placeholder.setPos(-50,180,5)
		self.sphereNodep.instanceTo(placeholder)

		self.taskMgr.add(self.update, "update")

	def update(self, task):
		global pos, vel
		#dt = globalClock.getDt()
		pos, vel = rk4(pos,vel,setdt)
		for i in range(len(self.spheres)):
			#print(pos[i])
			self.sphereNodep.getChild(i).setPos(pos[i,0],pos[i,1],pos[i,2])
			# dance in a sine pattern:
			#self.thetas[i] += (i/100 + 2.0) * dt
			#self.sphereNodep.getChild(i).setPos(i/4 - 50,180,(sin(self.thetas[i])*6)+10) 
			#self.sphereNodep.getChild(i).setHpr(cos(self.thetas[i]),0,sin(self.thetas[i]))

		return task.cont

app = TestBase()

#shader = Shader.load(Shader.SL_GLSL,
#                     vertex="testVertexShader.vert",
#                     fragment="testFragShader.frag")
#render.setShader(shader)

app.run()
