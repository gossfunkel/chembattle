from ursina import *
from ursina.prefabs.first_person_controller import *
import sys
import numpy as np
import atoms

WINDOW_WIDTH, WINDOW_HEIGHT = 1080, 600
app = Ursina(size=(WINDOW_WIDTH,WINDOW_HEIGHT))
window.title = "ChemBattle"
window.borderless = False
window.fullscreen = False
window.fps_counter.enabled = True
window.exit_button.enabled = False
window.color = color.black

p1ground = Entity(model='plane',
				scale=50,
				scale_z=100,
				origin_x=.5,
				origin_y=.5,
				texture='white_cube',
				texture_scale=(20,40),
				collider='mesh')
p2ground = Entity(model='plane',
				scale=50,
				scale_z=100,
				origin_x=-.5,
				origin_y=.5,
				texture='white_cube',
				texture_scale=(20,40),
				collider='mesh')

p1membrane = Entity(model='plane',
				scale=50,
				scale_z=100,
				rotation_z=-90,
				rotation_y=180,
				origin_y=.5,
				texture='white_cube',
				color=color.pink,
				alpha=0.2,
				double_sided=True,
				texture_scale=(20,10))
p2membrane = Entity(model='plane',
				scale=50,
				scale_z=100,
				origin_y=.5,
				rotation_z=90,
				rotation_y=180,
				texture='white_cube',
				color=color.yellow,
				alpha=0.2,
				double_sided=True,
				texture_scale=(20,10))

p1nucleus = Entity(origin=(4,2.5,3), 
					model='sphere',
					color=color.pink,
					collider='sphere',
					scale=10)
p2nucleus = Entity(origin=(-4,2.5,-3), 
					model='sphere',
					color=color.yellow,
					collider='sphere',
					scale=10)

camera.position = (0,150,0)
camera.rotation_x = 90
# this means
#           	 up = y-
# 				 /|\			   
# left - - - - - - - - - - - - - - right
# = x+			\|/					= x-
#           down = y+
# 
#	       _  front = z-
#  		   /| 	
#  back   /
#  = z+ |/_

def switchCameraPos():
	if camera.position != (0,150,0):
		camera.position = (0,150,0)
		camera.rotation_x = 90
		camera.rotation_y = 0
	else: 
		camera.position = (0,1.5,-175)
		camera.rotation_x = 0
		camera.rotation_y = 0

def input(self):
	if (held_keys['h']):
		print("adding Hydrogen")
		atoms.CreateAtom('Hydrogen')
	if (held_keys['j']):
		print("adding Carbon")
		atoms.CreateAtom('Carbon')
	if (held_keys['k']):
		print("adding Oxygen")
		atoms.CreateAtom('Oxygen')
	if (held_keys['l']):
		print("adding Nitrogen")
		atoms.CreateAtom('Nitrogen')

	if (held_keys['tab']):
		print("Switching views")
		switchCameraPos()

def update():
	camera.position += (5 * (held_keys['up_arrow'] - held_keys['down_arrow']) * time.dt, 0, 
						5 * (held_keys['up_arrow'] - held_keys['down_arrow']) * time.dt)
	camera.world_position += (0,0,20 * (held_keys['x'] - held_keys['z']) * time.dt) # z coord - z & x
	camera.world_position += (20 * (held_keys['d'] - held_keys['a']) * time.dt,0,0) # x coord - a & d
	camera.world_position += (0,20 * (held_keys['s'] - held_keys['w']) * time.dt,0) # y coord - w & s
	#camera.position_y += 20 * (held_keys['c'] - held_keys['v']) * time.dt
	camera.rotation_x += 20 * (held_keys['q'] - held_keys['e']) * time.dt # x rotation - q & e
	camera.rotation_y += 20 * (held_keys['v'] - held_keys['c']) * time.dt # y rotation - c & v

if __name__ == "__main__":
	app.run()

#cleanup
sys.exit()