from ursina import *
import sys
import numpy as np

WINDOW_WIDTH, WINDOW_HEIGHT = 1080, 600
app = Ursina(size=(WINDOW_WIDTH,WINDOW_HEIGHT))
window.title = "Card Battle"
window.borderless = False
window.fullscreen = False
window.fps_counter.enabled = True
window.exit_button.enabled = False
window.color = color.black

def dragme():
	testbutt.dragging = True

testbutt = Button(scale=1, color=color.green, origin=[0,0], label="testme", on_click=dragme, dragging=False)

def input(self):
	if testbutt.dragging:
		testbutt.origin = mouse.origin

# def update():
# 	pass

#if __name__ == "__main__":
#	sys.exit()
#else:
app.run()

#sys.exit()

