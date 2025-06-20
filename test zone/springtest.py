# testing grounds for constraint solving using rk4

# it's taken me a bit to get to how rk4 is used in games design
# when constraints are handled properly, f = ma is subbed to f = m*dv/dt, and from this we can get f*dt=m*dv
# this is why impulses pop out; f*dt = J . does this mean that J = p ??? p = m*v sooo does impulse map to momentum?

# when dealing with constraints, our solvers effectively apply smoothly interpolated impulses to correct velocities
# if done correctly, this will yield contact velocities and torques from friction etc
# it can also properly handle whole-body forces like magnetism better than the point-calculations that forces yield
# i've not finished fixing this test script to run again, but I'm hoping that if I'm actually getting anywhere with this,
# I'll get a constraint-solving function that spits out locations and acceleration/velocity for us to use for physics,
# which can call the force-calculating procedures

# i'm still wondering if there's utility to having the dljp, coul, bond functions output the force in the form of impulses,
# and if I can work those into our ODE then it might become implicit in the above substitution of Newton's 2nd law. 
# I think that gaining any sort of resolution from this might involve finding the appropriate way 

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

setdt = 0.02

block = Entity(model='cube',scale=3,color=color.blue, world_position=Vec3(0,10,0))
m = .1 		# mass in g
k = 1.0 		# spring constant
damp = 0.1 		# damping factor
f = np.zeros(3) # force
v = np.array([1.0,1.0,1.0]) # velocity 
t = 0.0 		# initial timestep
# appliedf = np.zeros(3) # force from input

# def input(self):
# 	global appliedf
# 	# press space to hit box
# 	if (held_keys['space']):
# 		appliedf = Vec3(0.0,100.0,0.0)
# 		window.color = color.brown
# 	else: 
# 		appliedf = Vec3(0.0,0.0,0.0)
# 		window.color = color.gray

# Applies more force the further the mass is from the equilibrium point k**2
#def spring(t,pos):
#	# dv/dt = a
#	#if t > 0:
#		return np.array([0,pos[1]*-(k**2) / dt,0]) MATHS CHECK
#	#else: 
#	#	return np.zeros(3)

def spring(pos,v): return (-k * pos) / damp * v

# gradient * change over x axis + integral constant
# TODO swap these out for my ODEs and correct the functioning of the system
# this doesnt do angular momentum or anything but i found a good page on my phone that has code i just need to translate from C(++?)
# maybe it should return a velocity (change in position) rather than a force?
def accel(dt,r): return (-k * r) * dt
def velacc(dt,a): return a / dt
def posvel(dt,v): return v / dt

# integrates by using values from the differential function dy/dx to estimate more accurate values for y
def rk4(x,y,dx,evaluate):
	# runge-kutte 4
	k1 = dx * evaluate(x,y) 				# y(x) at x
	k2 = dx * evaluate(x+dx/2.0,y+k1/2.0) 	# y(x) at x+dx/2
	k3 = dx * evaluate(x+dx/2.0,y+k2/2.0) 	# y(x) at x+dx/2
	k4 = dx * evaluate(x+dx,y+k3)			# y(x) at x+dx
	y = y + (k1 + 2*k2 + 2*k3 + k4)/6 		# weighted average for y(x) at x+dx
	x = x + dx
	return x,y

def update():
	global v
	global appliedf
	global block
	global t
	
	a = np.zeros(3) # accelerationa = np.zeros(3) # acceleration
	
	#appliedf += rk4step(time.dt,k,block.world_position,spring)			# spring
	#a = (-k* Vec3(0,0,0)-block.world_position)/ m
	# using spring constant to define a, calculate v using rk4 integration
	v = rk4(t,v,time.dt,spring)[1]
	#v = rk4(time.dt,v,time.dt,velacc)[1]
	# how to work out r when v(dt/2) is unknown? can i get rk4 to do both calculations at once?
	t, block.world_position = rk4(t,v,time.dt,posvel)
	#v = v + a/setdt 						# move
	#block.world_position = block.world_position + v/setdt

	print(block.world_position)

	# increment time
	#t += time.dt

if __name__ == "__main__":
	app.run()

sys.exit()
