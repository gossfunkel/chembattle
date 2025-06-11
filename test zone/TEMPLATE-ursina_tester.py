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
window.color 				= color.black

def update():
	pass

if __name__ == "__main__":
	app.run()

sys.exit()
