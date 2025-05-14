from ursina import *
import sys
import atoms
import AdvancedMoleculePhysics as amp
import numpy as np

#font = pygame.font.SysFont('Arial', 22)
WINDOW_WIDTH, WINDOW_HEIGHT = 1080, 600
menuOpenState = False
#textSurfs   = pygame.sprite.LayeredDirty()
#debugOutput = pygame.sprite.LayeredDirty()
uiobjects   = []
menubuttons = []
#all_sprites = pygame.sprite.Group()
Button.default_color = '#aaaaaa'

def AddUIObj(object):
	uiobjects.append(object)
	#all_sprites.add(object)

#def AddSprite(spr):
#	all_sprites.add(spr)

def GetUIObjects():
	return uiobjects

#def GetAllSprites():
#	return all_sprites

class CameraManager():
	def __init__(self):
		self.cameraFollowing = False

	def toggleFollow(self):
		self.cameraFollowing = not self.cameraFollowing
		if self.cameraFollowing:
			print("camera is following atom")

	def isFollowing(self):
		return self.cameraFollowing

cm = CameraManager()

def toggleMenu():
	# expand menu if closed, close menu if open
	#self.visible = (menuOpenState or not self.menuItem)
	#print(menuOpenState)
	global menuOpenState
	menuOpenState = not menuOpenState
	for butt in menubuttons:
		butt.disabled = not menuOpenState
		butt.visible = menuOpenState
	#if menuOpenState: menubuttons.draw(screen) #(menuOpenState or not self.menuItem)
	#else: menubuttons.clear(screen,'#020015')

def quitButt():
	# close game when clicked, if menu is open
	global menuOpenState
	if menuOpenState: 
		quit()
		sys.exit()

def createMol(atom):
	amp.CreateMolecule([0,0,0], atom.velocity, atom)

def createAt():
	global spawn
	atoms.CreateAtom(spawn)

class CustomButton(Button):
	def __init__(self, position, label, disabled=False, menuItem=False):
		super().__init__(parent=camera.ui,
						 text=label,
						 origin=position,
						 billboard=True,
						 highlight_color=self.color.tint(.1),
						 pressed_color=self.color.tint(-.2),
						 scale=0.1)
		self.menuItem = menuItem
		#print(self)
		#AddUIObj(self)
		if menuItem: 
			menubuttons.append(self)

class FollowButton(CustomButton):
	def __init__(self):
		super().__init__((-6, 4.2, 0), "FOLLOW")
		self.on_click = cm.toggleFollow
	
class MenuButton(CustomButton):
	def __init__(self):
		super().__init__((-8, 4.2, 0), "MENU")
		self.on_click = toggleMenu
		self.tooltip = Tooltip("Open the menu")

class QuitButton(CustomButton):
	def __init__(self):
		super().__init__((-8, 3, 0), "QUIT", menuItem=True)
		self.visible = False
		self.disabled = True
		self.on_click = quitButt
		self.tooltip = Tooltip("CAUTION: SAVING NOT IMPLEMENTED. ALL WILL BE LOST")

class AtomButton(CustomButton):
	def __init__(self):
		super().__init__((self.xco, self.yco, 0), ("Add\n" + self.nam))
		self.on_click = createAt

class HydrogenButton(AtomButton):
	def __init__(self):
		#print("hydrobutt initialising")
		self.xco = 7.8
		self.yco = 4.2
		self.nam = "Hydrogen"
		super().__init__()
		self.tooltip = Tooltip("Spawn a Hydrogen atom")
		global spawn 
		spawn = 'Hydrogen'


class CarbonButton(AtomButton):
	def __init__(self):
		self.xco = 6.5
		self.yco = 4.2
		self.nam   = "Carbon"
		super().__init__()
		self.tooltip = Tooltip("Spawn a Carbon atom")

class NitrogenButton(AtomButton):
	def __init__(self):
		self.xco = 5.2
		self.yco = 4.2
		self.nam   = "Nitrogen"
		super().__init__()
		self.tooltip = Tooltip("Spawn a Nitrogen atom")

class OxygenButton(AtomButton):
	def __init__(self):
		self.xco = 3.9
		self.yco = 4.2
		self.nam   = "Oxygen"
		super().__init__()
		self.tooltip = Tooltip("Spawn an Oxygen atom")
