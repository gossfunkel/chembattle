from direct.showbase.ShowBase import ShowBase
from panda3d.core import Shader, Geom, GeomNode, GeomTriangles, GeomVertexWriter
from panda3d.core import GeomVertexFormat, GeomVertexData, TransparencyAttrib
from math import sin, cos

def createColouredRect(x, z, width, height, colors=None):
	_format = GeomVertexFormat.getV3c4()
	vdata   = GeomVertexData('square', _format, Geom.UHDynamic)

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

class TestBase(ShowBase):
	def __init__(self):
		ShowBase.__init__(self)
		# get nodepath for environment model, then reparent to render node
		#self.scene = self.loader.loadModel("models/environment")
		#self.scene.reparentTo(self.render)
		# Apply scale and position transforms on the model.
		#self.scene.setScale(0.25)
		#self.scene.setPos(0, 42, 0)

		self.cam.setPos(-5, 0, 3)

		newSquare = createColouredRect(0, 0, 100, 100)

		geomNode = GeomNode('square')
		geomNode.addGeom(newSquare)
		self.render.attachNewNode(geomNode)

		# get nodepath for box, then reparent to render node
		self.sphere = self.loader.loadModel("models/misc/sphere")
		self.sphere.setScale(0.45)
		self.sphere.reparentTo(self.render)
		self.theta = 0.0

		self.taskMgr.add(self.update, "update")

	def update(self, task):
		dt = globalClock.getDt()
		self.theta += 2.0 * dt
		self.sphere.setPos((cos(self.theta)*8)-8,(sin(self.theta)*8)-8,2)
		self.sphere.setHpr(cos(self.theta),0,sin(self.theta))

		return task.cont

app = TestBase()


shader = Shader.load(Shader.SL_GLSL,
                     vertex="testVertexShader.vert",
                     fragment="testFragShader.frag")
camera.setShader(shader)

app.run()
