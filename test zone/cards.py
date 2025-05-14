import atoms
import MoleculePhysics as mp

# TODO : for all of these functions:
# find a way to handle the cards independently of the other player
# i think this will be sorted when i do the networking implementation

def Drag(self, loc, player):
	# drag and move the card around
	pass

def PlayerMoveDirectToggle(player, dirplayerDirect=0):
	#switch player direction of motion
	if playerDirect == 0:
		playerDirect=1
	else:
		playerDirect=0

def F_GrowFlaggelum(self, player):
	# Grow a flaggelum that allows you to select to move further from or closer to the other player by 1 map unit
	player.mobile=True
	self.on_click = (PlayerMoveDirectToggle, args=player)

class Flaggelum(Card):
	def __init__(self, player):
		super().__init__()
		self.activated = F_GrowFlaggelum

class Card(Button):
	def __init__(self, player):
		super().__init__(model='plane',
						 scale=2,
						 billboard=True)
		self.player = player
		self.on_click = drag

def F_DropofWater():
	# generate 50 water molecules in the appropriate player's zone
	mp.createMols('H20', np.zeros(3), 50, player)

class DropofWater(Card):
	def __init__(self, player):
		super().__init__(player)
		self.activated = F_DropofWater
		self.playToBoard = False;

def F_AcidSplash():
	# generate 15 simple acid molecules in the appropriate player's zone
	mp.createMols('SO4', np.zeros(3), 15, player)

class AcidSplash(Card):
	def __init__(self, player):
		super().__init__(player)
		self.activated = F_AcidSplash
		self.playToBoard = False;

def F_DropofLye():
	# generate 15 simple basic molecules in the appropriate player's zone
	mp.createMols('NaOH', np.zeros(3), 15, player)

class DropofLye(Card):
	def __init__(self, player):
		super().__init__(player)
		self.activated = F_DropofLye
		self.playToBoard = False;

def F_ElectronDeathRay():
	# generate a chlorophyll
	mp.createMols('chlorophyll', np.zeros(3), 1, player)

class ElectronDeathRay(Card):
	def __init__(self, player):
		super().__init__()
		self.activated = F_ElectronDeathRay
		self.playToBoard = False;

def F_PrimordialSoup():
	# give a handful of nucleotides and amino acids 
	mp.createMols('guanine', 0, 4, player)
	mp.createMols('adenine', 0, 4, player)
	mp.createMols('cytosine', 0, 4, player)
	mp.createMols('thymine', 0, 4, player)
	mp.createMols('pteridine', 0, 2, player)
	mp.createMols('purine', 0, 2, player)
	mp.createMols('uracil', 0, 4, player)
	mp.createMols('guanine', 0, 4, player)

class PrimordialSoup(Card):
	def __init__(self, player):
		super().__init__()
		self.activated = F_PrimordialSoup
		self.playToBoard = False;

def F_FirinMahLazah():
	uvbeam.on(3)

class FirinMahLazah(Card):
	def __init__(self, player):
		super().__init__()
		self.activated = F_FirinMahLazah
		self.playToBoard = False;

def F_DryingSalts():
	# provide some dehydrating salts
	mp.createMols('MgCO3',0,5,player)
	mp.createMols('CaCO3',0,5,player)
	mp.createMols('MgSO4',0,5,player)

class DryingSalts(Card):
	def __init__(self, player):
		super().__init__()
		self.activated = F_DryingSalts
		self.playToBoard = False;

def F_Glycolysis(inmol):
	# converts Glucose to Pyruvate
	if (inmol.name = 'Glucose'):
		mp.createMols('H2O',self.spawnpos,2,player)
		mp.createMols('H+',spawnpos,2,player)
		return mp.enzyme([inmol,'NAD+','NAD+','ADP','ADP'],['Pyruvate','NADH','NADH','ATP','ATP'])

class Glycolysis(Card):
	def __init__(self, player):
		super().__init__()
		self.enzyme = F_Glycolysis
		self.playToBoard = True;

def F_Lypogenesis(inmols):
	# converts Glucose to Pyruvate
	# TODO CHECK IF INPUTS ARE CORRECT
	return mp.enzyme([inmols[0],inmols[0],inmols[0]],['*Triglyceride'])

class Lypogenesis(Card):
	def __init__(self, player):
		super().__init__()
		self.enzyme = F_Lypogenesis
		self.playToBoard = True;
