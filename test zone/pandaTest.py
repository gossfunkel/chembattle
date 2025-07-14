from direct.showbase.ShowBase import ShowBase
from panda3d.core import Shader, Geom, GeomNode, GeomTriangles, GeomVertexWriter, loadPrcFileData
from panda3d.core import GeomVertexFormat, GeomVertexData, TransparencyAttrib, NodePath
from math import sin, cos
import numpy as np

def createColouredRect(x, z, width, height, colors=None):
	format = GeomVertexFormat.getV3c4()
	vdata   = GeomVertexData('square', format, Geom.UHDynamic)

	vertex  = GeomVertexWriter(vdata, 'vertex')
	color   = GeomVertexWriter(vdata, 'color')

	vertex.addData3(x, 0, z)
	vertex.addData3(x+width, 0, z)
	vertex.addData3(x+width, 0, z+height)
	vertex.addData3(x, 0, z+height)

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
  
loadPrcFileData('', 'win-size 1000 600') 
loadPrcFileData('', 'show-frame-rate-meter 1')

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
		sphereNum = 500 
		self.spheres = np.zeros(sphereNum)
		self.thetas = np.zeros(sphereNum)
		self.sphere = self.loader.loadModel("../sphere")
		self.sphereNodep = NodePath('spheres') # placeholders can be attached to another node to spawn groups of instances
		for i in range(sphereNum):
			placeholder = self.sphereNodep.attachNewNode("Sphere-Placeholder")
			placeholder.setScale(0.5)
			placeholder.setPos(0,i/4,0)
			self.sphere.instanceTo(placeholder)
			self.thetas[i] = 0.0
		placeholder = render.attachNewNode("SphereGroup")
		placeholder.setPos(-50,180,5)
		self.sphereNodep.instanceTo(placeholder)

		self.taskMgr.add(self.update, "update")

	def update(self, task):
		dt = globalClock.getDt()
		for i in range(len(self.spheres)):
			self.thetas[i] += (i/100 + 2.0) * dt
			self.sphereNodep.getChild(i).setPos(i/4 - 50,180,(sin(self.thetas[i])*6)+10)
			#self.sphereNodep.getChild(i).setHpr(cos(self.thetas[i]),0,sin(self.thetas[i]))

		return task.cont

app = TestBase()

#shader = Shader.load(Shader.SL_GLSL,
#                     vertex="testVertexShader.vert",
#                     fragment="testFragShader.frag")
#render.setShader(shader)

app.run()
