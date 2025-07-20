import sys
from direct.showbase.ShowBase import ShowBase
from panda3d.core import Shader, Geom, GeomNode, GeomTriangles, GeomVertexWriter, loadPrcFileData
from panda3d.core import GeomVertexFormat, GeomVertexData, TransparencyAttrib
from panda3d.core import NodePath, DirectionalLight, PointLight, LightRampAttrib
from direct.filter.CommonFilters import CommonFilters
from Ursina import Entity
#from math import sin, cos
import numpy as np
import PhysicsEngine as phys
import atoms
import ui
  
loadPrcFileData('', 'win-size 1080 600') 
loadPrcFileData('', 'show-frame-rate-meter 1')
loadPrcFileData('', 'hardware-animated-vertices true')
loadPrcFileData('', 'basic-shaders-only false')
loadPrcFileData('', 'window-title ChemBattle')
#loadPrcFileData('', '')


p1ground = Entity(model='plane',
				world_position=(25,0,25),
				scale_x=50,
				scale_z=50,
				texture='white_cube',
				texture_scale=(50,40),
				collider='mesh')
p2ground = Entity(model='plane',
				world_position=(25,0,75),
				scale_x=50,
				scale_z=50,
				texture='white_cube',
				texture_scale=(50,40),
				collider='mesh')

p1membrane = Entity(model='plane',
				world_position=(25,5,20),
				scale=50,
				scale_y=10,
				scale_z=50,
				rotation_z=90,
				rotation_y=90,
				texture='white_cube',
				color=color.pink,
				alpha=0.15,
				double_sided=True,
				collider='box',
				texture_scale=(20,10))
p2membrane = Entity(model='plane',
				world_position=(25,5,80),
				scale=50,
				scale_y=10,
				scale_z=50,
				#origin_y=7.5,
				rotation_z=90,
				rotation_y=90,
				texture='white_cube',
				color=color.yellow,
				alpha=0.15,
				double_sided=True,
				collider='box',
				texture_scale=(20,10))

p1nucleus = Entity(world_position=(10,0,8), 
					model='sphere',
					color=color.pink,
					collider='sphere',
					scale=10)
p2nucleus = Entity(world_position=(40,0,92), 
					model='sphere',
					color=color.yellow,
					collider='sphere',
					scale=10)

# INITIALISE THE UI
hydrobutt = ui.HydroButton()
carbobutt = ui.CarbondioxButton()
ammobutt  = ui.AmmoniaButton()
oxybutt   = ui.OxygenButton()
menbutt   = ui.MenuButton()
quitbutt  = ui.QuitButton()
follbutt  = ui.FollowButton()
#uiobjects.append()

#camera.rotation_y = -90
# this means
#           	 up = y+
# 				 /|\			   
# left - - - - - - - - - - - - - - right
# = x-			\|/					= x+
#           down = y-
# 
#	       _  back = z+
#  		   /| 	
#  front  /
#  = z- |/_

def switchCameraPos():
	if camera.position != (0,150,0):
		camera.position = (0,150,0)
		camera.rotation_x = 90
		camera.rotation_y = 0
	else: 
		camera.position = (25,50,-65)
		camera.rotation_x = 0
		camera.rotation_y = 0

def input(self):
	#if (held_keys['h']):
	#	atoms.CreateAtom('Hydrogen')
	#if (held_keys['j']):
	#	atoms.CreateAtom('Carbon')
	#if (held_keys['k']):
	#	atoms.CreateAtom('Oxygen')
	#if (held_keys['l']):
	#	atoms.CreateAtom('Nitrogen')

	if (held_keys['tab'] and not ui.cm.isFollowing()):
		print("Switching views")
		switchCameraPos()

class Chembattle(ShowBase):
	def __init__(self):
		ShowBase.__init__(self)

		# position camera
		self.cam.setPos(25,65,-70)
		self.cam.setHpr(30,0,0)

		# initialise shading and filters
		render.setShaderAuto()
		self.filters = CommonFilters(base.win,base.cam)
		self.filters.setBloom(blend=(0.25, 0.22, 0.28, 0.0), maxtrigger=1.1, desat=0.7, size='large')

		# add lights

		# instantiate ui objects

		self.taskMgr.add(self.update, "update", taskChain='default')

	def update(self, task):
		atomPos, elecPos = phys.update() 
		# which unit is handling rendering? PhysicsEngine.update() should probably return the positions

		if not ui.cm.isFollowing():
			camera.position += (5 * (held_keys['up_arrow'] - held_keys['down_arrow']) * time.dt, 0, 
								5 * (held_keys['up_arrow'] - held_keys['down_arrow']) * time.dt)
			camera.world_position += (0,0,20 * (held_keys['x'] - held_keys['z']) * time.dt) # z coord - z & x
			camera.world_position += (20 * (held_keys['d'] - held_keys['a']) * time.dt,0,0) # x coord - a & d
			camera.world_position += (0,20 * (held_keys['s'] - held_keys['w']) * time.dt,0) # y coord - w & s
			#camera.position_y += 20 * (held_keys['c'] - held_keys['v']) * time.dt
			camera.rotation_x += 20 * (held_keys['q'] - held_keys['e']) * time.dt # x rotation - q & e
			camera.rotation_y += 20 * (held_keys['v'] - held_keys['c']) * time.dt # y rotation - c & v
		else: 
			camera.look_at(mols[0].world_position)
			dist = np.linalg.norm(camera.world_position-mols[0].world_position)
			if dist > 10:
				amp.SlideTo(camera.world_position, mols[0].world_position, 15)

		#p1memcoll = p1membrane.intersects()
		#if p1memcoll.hit:

if __name__ == "__main__":
	app.run()

#cleanup
sys.exit()