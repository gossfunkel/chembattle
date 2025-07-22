from direct.showbase.ShowBase import ShowBase
from panda3d.core import Shader, Geom, GeomNode, GeomTriangles, GeomVertexWriter
from panda3d.core import GeomVertexFormat, GeomVertexData, TransparencyAttrib
from panda3d.core import NodePath, DirectionalLight, PointLight, LightRampAttrib
from panda3d.core import loadPrcFileData
from direct.interval.IntervalGlobal import *
from direct.interval.IntervalManager import ivalMgr
from direct.filter.CommonFilters import CommonFilters
from math import sin, cos
import numpy as np

#from panda3d.core import Thread

loadPrcFileData('', 'win-size 1400 900') 
loadPrcFileData('', 'show-frame-rate-meter 1')
loadPrcFileData('', 'hardware-animated-vertices true')
loadPrcFileData('', 'basic-shaders-only false')
loadPrcFileData('', 'threading-model Cull/Cull') # creates a different two-thread model: App is on its own thread, 
										# and Cull and Draw are together on a separate thread. This is most appropriate
										# when the total time doing App is similar to the total time for Cull + Draw.

	## how many spheres
sphereNum = 20
	## how many electrons
elecNum = 20

	### some physical constants
kb  = 0.8314459920816467 				  # Boltzmann constant
NA  = 6.0221409e+26 					  # Avogadro's constant x 1000 (g->kg)
ech = -1.60217662E-19 					  # electron charge (coulombs)
pch = 1.60217733E-19 					  # proton charge (coulombs)
kc  = 8.9875517923E9*NA*1E30*ech*ech/1E24 # electrostatic const (Daltons, electron charge, picosecond, angstrom units)
amH = 1.6735575E-27 					  # atomic mass of Hydrogen (kg)
gravConst = 6.6743e-11 					  # G (m3/kg/s)
	### physical values
chrg = np.array([(pch + ech) for _ in range(sphereNum)]) 	# TODO: ech moves with electrons
mass = np.array([amH for _ in range(sphereNum)])
#chrg = np.array([0.41 for _ in range(sphereNum)]) 			# what units are these? how could a hydrogen atom have 4100000000000000000 electrons' worth of charge???
#mass = np.array([0.001 for _ in range(sphereNum)]) 		# units ?? if these are kg, each 'hydrogen' atom weighs a gram
	## constant deltatime - n.b. for scientific accuracy, we require this to have a resolution of under 2fs (2e-15 seconds). 
		# The current value is: 0.1as! (1e-19 secs) 
		# Appropriately scaled masses and charges are required to get to realistic timescales
		# n.b. a ray of light travels the bond length of molecular hydrogen (H2) in 0.247attoseconds
		# atomic unit of time is 24.189as. neutral pions decay into gamma rays after 85as.
setdt = 0.1E-18
	## derivatives of motion for calculations
acc = np.zeros((sphereNum,3)) 									# a in m/s/s? A/ps/ps?
vel = np.zeros((sphereNum,3)) 									# v in m/s? A/ps?
#pos = np.array([[(sin(i)*2.)+50.,cos(i)*2.0-150.0,-15.0-i/4.] for i in range(sphereNum)])
pos = np.array([[50.+i,-150.,-15.] for i in range(sphereNum)])  # r in m? Angstroms?
electrons = np.array([pos]) 									# electron positions (m? A?)
#ep = np.array([3.24 for _ in range(sphereNum)])
#sg = np.array([0.98 for _ in range(sphereNum)])
	## Lennard-Jones parameters for Hydrogen's dipole effects
ep = 3.24
sg = 0.98

def createColouredRect(x, y, z, width, height, colors=None):
 	format = GeomVertexFormat.getV3n3c4()
 	vdata   = GeomVertexData('square', format, Geom.UHDynamic)

 	vertex  = GeomVertexWriter(vdata, 'vertex')
 	color   = GeomVertexWriter(vdata, 'color')
 	normal 	= GeomVertexWriter(vdata, 'normal')

 	vertex.addData3(x, y, z)
 	vertex.addData3(x+width, y, z)
 	vertex.addData3(x+width, y, z+height)
 	vertex.addData3(x, y, z+height)

 	normal.addData3(x*x, y*(y-1), z*z)
 	normal.addData3((x+width)*(x+width), y*(y-1), z*z)
 	normal.addData3((x+width)*(x+width), y*(y-1), (z+height)*(z+height))
 	normal.addData3(x*x, y*(y-1), (z+height)*(z+height))

 	if colors:
 		if len(colors) < 4:
 			colors = (colors[0],colors[0],colors[0],1.0)
 		color.addData4f(colors)
 		color.addData4f(colors)
 		color.addData4f(colors)
 		color.addData4f(colors)
 	else:
 		color.addData4f(1.0,0.0,0.0,0.0)
 		color.addData4f(0.0,1.0,0.0,0.0)
 		color.addData4f(0.0,0.0,1.0,0.0)
 		color.addData4f(0.0,0.0,0.0,1.0)

 	tris = GeomTriangles(Geom.UHDynamic)
 	tris.addVertices(0,1,2)
 	tris.addVertices(2,3,0)

 	square = Geom(vdata)
 	square.addPrimitive(tris)
 	return square

# dLJP for Hydrogen
def dLJP(dr,drv):
	#print(parti)
	dLJP=0.0
	# can I swap this for some kind of threshold function to avoid the branch?
	if ((dr[i] < 3) for i in dr): 												# calculate dLJP (8-14 from 6-12)

		r8vs = np.sum(np.transpose(np.transpose(drv) * ep*(sg**6) * (1.0/np.array(dr))**8),axis=0)
		r14vs= np.sum(np.transpose(np.transpose(drv) * 2.0*ep*(sg**12) * (1.0/np.array(dr))**14),axis=0)
		dLJP = (r14vs-r8vs)*24.0
	return dLJP

def coul(dr,drv,i,qi=0.41,qj=0.41):
	r3  = kc * qi * qj * ((1.0/np.array(dr))**3.0) 								# Coulomb's law BROKEN ONLY WORKS FOR TWO CHARGES
	return np.sum(np.transpose(np.transpose(drv)*r3), axis=0)					# transpose the distance vector, multiply by force, transpose back

def grav(dr,i):
	return np.sum(gravConst * mass[i] * mss / np.transpose(np.transpose(dr)**2) for mss in mass)

class TestBase(ShowBase):
	def __init__(self):
		ShowBase.__init__(self)

		render.setShaderAuto()

		self.filters = CommonFilters(base.win,base.cam)
		self.filters.setBloom(blend=(0.25, 0.22, 0.28, 0.0), maxtrigger=1.1, desat=0.7, size='large')

		self.cam.setPos(2, 25, -6)
		self.cam.setHpr(0,-40,0)

		dirLight 	= DirectionalLight('dirLight')
		dirLight.setColorTemperature(6650)
		dirLight.setShadowCaster(True, 512, 512)
		dirLightNp  = render.attachNewNode(dirLight)
		dirLightNp.setHpr(45,-50,0)
		render.setLight(dirLightNp)

		pointLight 	= PointLight('poiLight')
		pointLight.setColor((0.2, 0.6, 1.1, 1))
		pointLight.setShadowCaster(True, 512, 512)
		#pointLight.attenuation = (1,0,1)
		plnp = render.attachNewNode(pointLight)
		plnp.setPos(40, 100, 120)
		render.setLight(plnp)

		render.setAttrib(LightRampAttrib.makeHdr0())

		newSquare = createColouredRect(-150, 150, -250, 300, 250)

		geomNode = GeomNode('square')
		geomNode.addGeom(newSquare)
		self.render.attachNewNode(geomNode)
		#geomNode.getParent(0).setLightOff()

		## this should now use hardware instancing to generate many copies of the one sphere model:
		self.spheres = np.zeros(sphereNum)
		self.elecs = np.zeros(elecNum)
		self.sphere = self.loader.loadModel("../sphere")
		self.sphere.setScale(0.5)
		self.sphereNodep = NodePath('spheres') 
		# placeholders can be attached to another node to spawn groups of instances
		for i in range(sphereNum):
			placeholder = self.sphereNodep.attachNewNode("Sphere-Placeholder")
			placeholder.setScale(0.5)
			placeholder.setPos(pos[i,0],pos[i,1],pos[i,2])
			self.sphere.instanceTo(placeholder)
			#self.thetas[i] = 0.0
		placeholder = render.attachNewNode("SphereGroup")
		placeholder.setPos(-50,180,5)
		self.sphereNodep.instanceTo(placeholder)

		# tasks method commented out. To reactivate, ensure to add 'async' to update method declaration
		#self.asyncTaskMgr = AsyncTaskManager("asyncTaskMgr").getGlobalPtr()
		# create a task chain capable of multithreading
		#self.taskMgr.setupTaskChain('physTaskChain', 
		#							numThreads=8, 
		#							threadPriority=1)
		# initialise variables and weights for rk4 in scope, each as an array of size sphereNum, populated with 3d vectors
		self.k = np.array([np.zeros((sphereNum,3)) for _ in range(4)])
		self.r = pos
		self.v = vel
		self.t = 0.0
		self.rk4step = 0
		#self.tasks = []
		# generate an interval to calculate for each sphere
		self.odeParallel = Parallel(name="odeParallel")
		for i in range(sphereNum):
			self.odeParallel.append(Func(self.ode, i))
			#self.tasks.append(taskMgr.doMethodLater(0, self.ode, "sphere_ode", extraArgs=[i], taskChain="physTaskChain", 
			#													 sort=0, appendTask=True))
		self.taskMgr.add(self.update, "update", taskChain='default')
	
	def ode(self, i): 
		drv = np.delete(self.r - self.r[i],i,0) 							# distance vectors from r[i] to all *other* positions
		dr  = [np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) for a in drv] 		# absolute distance of that lad
		force = -dLJP(dr,drv) + coul(dr,drv,i) #- grav(dr,i) 				# calculate F   - 	!!! CURRENTLY, ADDING GRAVITY HALVES THE FRAMERATE
		self.k[self.rk4step,i] = np.array([np.transpose(np.transpose(force) / mass[i])]) 	#!!!!!!!!!! hacky - only works while masses are equal
		self.v[i] = self.v[i] + self.k[self.rk4step,i] * self.t 			# find r'(t) = v(t) from a = r''(t)
		self.r[i] = self.r[i] + self.v[i] * self.t 							# find r(t)
		return 1 # task.done	

	def update(self, task):
		global pos, vel
		
		#dt = globalClock.getDt()
		#self.plnp.setPos(180*sin(dt), 100*sin(dt), 200)

		# run 4 physics calculations per tick
		for j in range(4):
			for i in range(4):
				#print("stepping into rk4 stage " + str(i+1))
				#currentTsk = self.taskMgr.getCurrentTask()
				#Task.sequence(taskMgr.getTasks())
				#if (taskMgr.mgr.getActiveTasks().hasTask(self.tasks[sphereNum-1])): # wait for previous round of tasks to finish
				#	for tsk in taskMgr.mgr.getActiveTasks():
				#		if (tsk.name == "sphere_ode"):
				#			#print("awaiting " + str(tsk))
				#			await tsk
				self.r = pos 
				self.rk4step = i
				if (i == 0): 
					self.t = 0.
					self.v = vel
				elif (i == 1 or i == 2): 
					self.t = setdt/2. 
					self.v = vel + self.k[i-1] * setdt/2
				else: 
					self.t = setdt
					self.v = vel + self.k[3] * setdt

				self.odeParallel.start()

				#for tsk in self.taskMgr.getDoLaters(): # run one task per sphere
				#	if (tsk.name == "sphere_ode"):
				#		#print(tsk)
				#		tsk.again
			# (after rk4, inside physics iterations)
			#if (taskMgr.getTasksNamed("sphere_ode") != 0): # wait for previous round of tasks to finish
			#	for tsk in taskMgr.getTasksNamed("sphere_ode"):
			#		#print("awaiting " + str(tsk))
			#		await tsk
			vel = self.v + setdt/6*(self.k[0] + 2*self.k[1] + 2*self.k[2] + self.k[3])
			pos = self.r
		#print(pos)
		
		# # update the positions of the spheres:
		#self.sphereNodep.getChildren().setPos(pos)
		for i in range(sphereNum):
			#print(pos[i])
			self.sphereNodep.getChild(i).setPos(pos[i,0],pos[i,1],pos[i,2])

		return task.cont

app = TestBase()

#shader = Shader.load(Shader.SL_GLSL,
#                     vertex="testVertexShader.vert",
#                     fragment="testFragShader.frag")
#render.setShader(shader)

app.run()
