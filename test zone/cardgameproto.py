# prototyping the card system, collectibles, and other ui/interactive/player-side features

from ursina import *
from ursina.prefabs.first_person_controller import *
from ursina.shaders import lit_with_shadows_shader 
from queue import Queue
#from typing import List
from random import uniform
from dataclasses import dataclass
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
EditorCamera()
pivot = Entity()
DirectionalLight(parent=pivot, y=2, z=3, shadows=True, rotation=(45, -55, 45))

atpStock = []
atpCount = 0
atpIcon = Button(color=color.cyan,scale=0.1,text_color=color.black,text='ATP:',position=(.8,.4),parent=camera.ui)
atpText = Text('0',position=(.8,.3735),text_color=color.black,background=color.gray,parent=camera.ui)

nucStock = []
nucCount = 0
nucIcon = Button(color=color.magenta,scale=0.1,text_color=color.black,text='NP:',position=(.8,.2),parent=camera.ui,radius=.5)
nucText = Text('0',position=(.8,.1735),text_color=color.black,background=color.gray,parent=camera.ui)

cardsprite = load_texture('../assets/card-tex-paper.jpg')

# smoothly transition a vector via linear interpolation
def SlideTo(pops, fpos, lerp):
	#target = Vec3([pops])
	#follow = Vec3([fpos])
	distance = np.linalg.norm(pops-fpos)
	if (distance > 0):
		direction = (pops - fpos) / distance
		min_step = max(0.02, distance - 0) #final int is max follow distance setting
		max_step = distance
		step_distance = min_step + (max_step - min_step) * lerp #LERP factor
		new_follow = Vec3(fpos + direction * (step_distance * time.dt))
	else: return(fpos)
	return(new_follow)

class Card(Draggable):
	def __init__(self, name, maxtemp, text, funcFx):
		super().__init__(text=(name + "\n\n\n" + text + "\n\n\nMax temp: " + str(maxtemp)),
			color=color.rgb(0.5,0.5,0),
			scale_x=0.3,
			scale_y=0.4,
			position=(.8,-0.3),
			texture=cardsprite,
			disabled=True,
			shader=lit_with_shadows_shader,
			parent=camera.ui,
			highlight_scale=1.05,
			radius=.04)
		self.name = name
		self.maxtemp = maxtemp
		self.fx = Func(funcFx)
		# lerp variables - should probably have a method to set lerpspeed for a movement rather than randomly changing it
		self.destination = self.position
		self.lerpspeed = 14.6

	def __repr__(self):
		return (f'''{self.__class__.__name__} named {self.name} with effect {self.fx}.''')

	def stop_dragging(self):
		super().stop_dragging()
		print("dropping card")
		if (mouse.position[0] > -.5 and mouse.position[0] < .5 and mouse.position[1] > -.5 and mouse.position[1] < .5):
			self.fx()
			discardCard(self)
		else:
			print(mouse.position)
			# todo change for reliable hand positions . some kind of hand tracker
			self.destination = Vec3(-.7 + counter/5,-.5,0)
			self.lerpspeed = 16.8

	def update(self):
		super().update()
		distancetodest = np.linalg.norm(self.destination - self.position)
		if ((not self.dragging) and distancetodest > 0.05): 
			print("lerping card")
			self.position = SlideTo(self.destination,self.position,self.lerpspeed)
		elif ((not self.dragging) and distancetodest < 0.05 and distancetodest > 0): 
			self.position = self.destination

class ATP(Entity):
	def __init__(self):
		super().__init__(self,
			model='cube',
			parent=scene,
			scale=0.2,
			collider='cube',
			color=color.cyan,
			origin=(-1 + uniform(0,1),1 + uniform(0,1)))
		self.destination = self.position
		self.lerpspeed = 20.0

	def update(self):
		global atpCount

		self.rotation_y += 10 * time.dt

		if self.hovered:
			self.destination = Vec3(6.5,3.2,0)
		distancetodest = np.linalg.norm(self.destination - self.position)
		if (distancetodest > 0.05): 
			print("lerping ATP")
			self.position = SlideTo(self.destination,self.position,self.lerpspeed)
		elif (distancetodest < 0.05 and distancetodest > 0): 
			self.enabled = False
			atpCount += 1

class Nucleotide(Entity):
	def __init__(self):
		super().__init__(self,
			model='sphere',
			parent=scene,
			scale=0.24,
			collider='sphere',
			color=color.magenta,
			origin=(1 + uniform(0,1),-1 + uniform(0,1)))
		self.destination = self.position
		self.lerpspeed = 20.0

	def update(self):
		global nucCount

		self.rotation_y += 10 * time.dt

		if self.hovered:
			self.destination = Vec3(6.5,1.8,0)
		distancetodest = np.linalg.norm(self.destination - self.position)
		if (distancetodest > 0.05): 
			print("lerping Nucleotide")
			self.position = SlideTo(self.destination,self.position,self.lerpspeed)
		elif (distancetodest < 0.05 and distancetodest > 0): 
			self.enabled = False
			nucCount += 1

def FXmakeATP(): return atpStock.append(ATP())
def FXmakeNuc(): return nucStock.append(Nucleotide())

@dataclass
class CardDeck():
	cards: list[Card]
	size: int = 0

	def __post_init__(self):
		self.size = len(self.cards)

testcardlist = []
for c in range(60):
	if c%2:
		testcardlist.append(Card('Free nucleotides',
				250.0,
				"Generate some Nucleo-\ntides in your virtual\nenvironment.\nOut of nothing.\nPoof.",
				Func(FXmakeNuc)))
	else: 
		testcardlist.append(Card('Free ATP',
				250.0,
				"Generate some ATP\nin your virtual\nenvironment.\nOut of nothing.\nPoof.",
				Func(FXmakeATP)))

testDeck = CardDeck(testcardlist)


activeDeck = Queue(maxsize=60)
cardHand = [] # hand of cards
counter = 0 

def loadDeck(deck):
	global activeDeck
	#print("loading " + str(deck))
	# check that deck has 60 cards in it
	if deck.size > 60:
		raise Exception("Deck too large! Must be exactly 60 cards.")
	if deck.size < 60:
		raise Exception("Deck too small! Must be exactly 60 cards.")

	# put each card from the deck in the active deck
	#deck.shuffle()
	for card in deck.cards:
		#print ("loading card " + str(card))
		activeDeck.put(card)
	print("deck finished loading.")
	#print(len(activeDeck))

def drawCard():
	global counter
	global cardHand
	if len(cardHand) < 10:
		counter += 1
		drawnCard = activeDeck.get()
		cardHand.append(drawnCard)
		drawnCard.destination = Vec3(-.7 + counter/5,-.5,0)
		drawnCard.lerpspeed = 19.0
		drawnCard.color=color.rgb(0.5 + counter/25,0.5 + counter/25,0 + counter/25)
		drawnCard.enabled = True
		return drawnCard
	else:
		print("Your hand is full of cards!")
		return False

def discardCard(disCard):
	global counter
	global cardHand
	if type(disCard) == None:
		# TODO throw error
		print("Can't discard; no card selected to discard!")
		return None

	# TODO double check the card is in the hand, throw error otherwise
	cardHand.remove(disCard)
	counter -= 1
	disCard.enabled = False

	# if not buttq.empty():
	# 	counter -= 1
	# 	buttq.get().enabled = False # does this actually remove the object from memory?
	# else:
	# 	print("Queue already empty!")
	# 	return False

drawbutt = Button(scale_x=.4, 
				scale_y=.12, 
				color=color.green, 
				position=(-.75,.4), 
				parent=camera.ui,
				text="Draw a card", 
				on_click=drawCard)
#discardbutt = Button(scale_x=.4, scale_y=.12, color=color.red, position=(-.75,.24), text="despawn a card", on_click=discardCard)

# temporary
loadDeck(testDeck)

#def input(self):
#	if 

def update():
	atpText.text = atpCount
	nucText.text = nucCount
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

