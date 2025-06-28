from direct.showbase.ShowBase import ShowBase
from panda3d.core import Shader#
from math import sin, cos

class TestBase(ShowBase):
	def __init__(self):
		ShowBase.__init__(self)
		# get nodepath for environment model, then reparent to render node
		self.scene = self.loader.loadModel("models/environment")
		self.scene.reparentTo(self.render)
		# Apply scale and position transforms on the model.
		self.scene.setScale(0.25)
		self.scene.setPos(0, 42, 0)

		self.cam.setPos(-5, 0, 3)

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
app.scene.setShader(shader)

app.run()
