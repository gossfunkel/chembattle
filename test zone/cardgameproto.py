from ursina import *
from queue import Queue
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

buttq = Queue(maxsize=10) # hand of cards
counter = 0 

def spawn():
	global counter
	if not buttq.full():
		counter += 1
		return buttq.put(Draggable(scale_x=0.3, scale_y=0.4, color=color.yellow, position=(-.7 + counter/5,-.5), text="testme"))
	else:
		print("Queue full!")
		return False

def despawn():
	global counter
	if not buttq.empty():
		counter -= 1
		buttq.get().enabled = False # does this actually remove the object from memory?
	else:
		print("Queue already empty!")
		return False

spawnbutt = Button(scale_x=.4, scale_y=.12, color=color.green, position=(-.75,.4), text="spawn a card", on_click=spawn)
rmvbutt = Button(scale_x=.4, scale_y=.12, color=color.red, position=(-.75,.24), text="despawn a card", on_click=despawn)

#def input(self):
#	pass

def update():
	pass
	#if buttq.full:
	#	despawn()
	#	window.color = color.orange
	#else:
	#	window.color = color.black

#if __name__ == "__main__":
#	sys.exit()
#else:
app.run()

sys.exit()

