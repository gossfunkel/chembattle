import numpy as np 
from ursina import *
from ursina.shaders import lit_with_shadows_shader
from sys import exit
import math


WINDOW_WIDTH, WINDOW_HEIGHT = 600, 600
app 						= Ursina(size=(WINDOW_WIDTH,WINDOW_HEIGHT))
EditorCamera()
pivot = Entity()
camera.fov = 100
DirectionalLight(parent=pivot,y=2,z=3,shadows=True,rotation=(45,-45,45))
window.borderless 			= False
window.fullscreen 			= False
window.fps_counter.enabled  = True
window.exit_button.enabled  = False
window.color 				= color.gray

setdt = 0.0002

block = Entity(model='cube',scale=3,color=color.blue, world_position=np.zeros(3))
m = 10.0 		# mass in g
k = 0.25 		# spring constant
f = np.zeros(3) # force
a = np.zeros(3) # acceleration
v = np.zeros(3) # velocity 
t = 0.0 		# initial timestep
appliedf = np.zeros(3) # force from input

def input(self):
	global appliedf
	# press space to hit box
	if (held_keys['space']):
		appliedf = Vec3(0.0,1.0,0.0)
		window.color = color.brown
	else: 
		appliedf = Vec3(0.0,0.0,0.0)
		window.color = color.gray


def secondLaw(t,v):
	# dv/dt = a
	if t > 0:
		return (f/m) / t # a = f/m and v = a/t
	else: 
		return 0

def spring(t,pos):
	# dr/dt = v
	if t > 0:
		return np.array([0,pos[1]*k**2 / t,0])
	else: 
		return np.zeros(3)

# def spring(t,vel):
# 	return np.array([vel[0],-k**2 * vel[1],vel[2]])

# def velacc(t,vel,a):
# 	vel = vel +a * t
# 	return vel

# def posvel(t,r):
# 	return block.world_position + v * t

def rk4step(t,dt,x,evaluate):
	# runge-kutte 4
	# f(x,t) = dx/dt
	k1 = evaluate(t,x)
	k2 = evaluate(t + setdt/2,x+k1*setdt/2)
	k3 = evaluate(t + setdt/2,x+k2*setdt/2)
	k4 = evaluate(t + setdt,x + k3*setdt)
	return setdt * (k1 + 2*k2 + 2*k3 + k4)/6

def update():
	global v
	global f
	global block
	global t
	
	f = appliedf - rk4step(t,setdt,block.world_position,spring)			# spring
	#a = f / m
	v = rk4step(t,setdt,block.world_position,secondLaw)		# apply force
	#v = v + a/setdt 						# move
	block.world_position = block.world_position + v/setdt

	print(block.world_position)

	# increment time
	t += setdt

if __name__ == "__main__":
	app.run()

sys.exit()
