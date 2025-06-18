from ursina import *
from ursina.prefabs.first_person_controller import *
from queue import Queue
#from typing import List
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

boxes = []
cardsprite = load_texture('../assets/card-tex-paper.jpg')

class Card(Draggable):
	def __init__(self, name, maxtemp, text, funcFx):
		super().__init__(text=(name + "\n\n\n" + text + "\n\n\nMax temp: " + str(maxtemp)),
			color=color.rgb(0.5,0.5,0),
			scale_x=0.3,
			scale_y=0.4,
			position=(.8,.1),
			texture=cardsprite,
			disabled=True,
			highlight_scale=1.05,
			radius=.04)
		self.name = name
		self.maxtemp = maxtemp
		self.fx = Func(funcFx)

	def __repr__(self):
		return (f'''{self.__class__.__name__} named {self.name} with effect {self.fx}.''')

	def stop_dragging(self):
		if (mouse.position[0] > -.5 and mouse.position[0] < .5 and mouse.position[1] > -.5 and mouse.position[1] < .5):
			self.fx()
			discardCard(self)
		else:
			# todo: lerp to position
			print(mouse.position)
			self.position = (-.7 + counter/5,-.5)

def FXmakeBox(): return boxes.append(Entity(model='cube',parent=scene,scale=0.2,color=color.cyan, origin=(0,0)))

@dataclass
class CardDeck():
	cards: list[Card]
	size: int = 0

	def __post_init__(self):
		self.size = len(self.cards)

boxDeck = CardDeck([Card('Make Box',
				250.0,
				"Generate a three-\ndimensional box in your\nvirtual environment.\nOut of nothing.\nPoof.",
				Func(FXmakeBox)) for c in range(60)])


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
		drawnCard.position = (-.7 + counter/5,-.5)
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

drawbutt = Button(scale_x=.4, scale_y=.12, color=color.green, position=(-.75,.4), text="Draw a card", on_click=drawCard)
#discardbutt = Button(scale_x=.4, scale_y=.12, color=color.red, position=(-.75,.24), text="despawn a card", on_click=discardCard)

# temporary
loadDeck(boxDeck)

#def input(self):
#	if 

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

