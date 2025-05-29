from ursina import *
import sys
import atoms
import AdvancedMoleculePhysics as amp
import numpy as np
import functools

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

def createMol(*mols):
	# TODO: UPDATE TO PASS NECESSARY ARGS
	if len(mols) == 1:
		print ("creating " + str(mols[0]) + " molecule.")
		amp.CreateMolecule([0,0,0], np.zeros(3), mols[0])
	else:
		print("ui.createMol requires 1 argument passed! Tell it what molecule you want to create")
		return False

def createAt(atom):
	print ("Creating " + str(atom) + " atom.")
	createMol(atoms.CreateAtom(atom))

class CustomButton(Button):
	def __init__(self, position, label, arg=None, disabled=False, menuItem=False):
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

class SpawnButton(CustomButton):
	def __init__(self):
		super().__init__((self.xco, self.yco, 0), ("Add\n" + self.nam), self.nam)
		self.on_click = functools.partial(createMol, self.nam)
		#self.args = (self.nam,)

class HydroButton(SpawnButton):
	def __init__(self):
		#print("hydrobutt initialising")
		self.xco = 7.8
		self.yco = 4.2
		self.nam = 'H2O'
		super().__init__()
		self.tooltip = Tooltip("Spawn a water molecule")
		#self.args = ('Hydrogen',)

class CarbondioxButton(SpawnButton):
	def __init__(self):
		self.xco = 6.5
		self.yco = 4.2
		self.nam   = 'CO2'
		super().__init__()
		self.tooltip = Tooltip("Spawn a carbon dioxide molecule")
		#self.args = "H20"

class AmmoniaButton(SpawnButton):
	def __init__(self):
		self.xco = 5.2
		self.yco = 4.2
		self.nam   = 'NH3'
		super().__init__()
		self.tooltip = Tooltip("Spawn an ammonium molecule")
		#self.args = "Nitrogen"

class OxygenButton(SpawnButton):
	def __init__(self):
		self.xco = 3.9
		self.yco = 4.2
		self.nam   = 'O2'
		super().__init__()
		self.tooltip = Tooltip("Spawn a molecule of diatomic Oxygen")
		#self.args = "Oxygen"
