import atoms
import MoleculePhysics as mp

def Drag():
	# 
	pass

class Card(Button):
	def __init__(self, player):
		super().__init__(model='plane',
						 scale=2,
						 billboard=True)
		self.player = player
		self.on_click = drag

def F_DropofWater():
	# generate 50 water molecules in the appropriate player's zone
	# TODO : for all of these functions:
	# find a way to handle the cards independently of the other player
	# i think this will be sorted when i do the networking implementation
	mp.createMols('H20', np.zeros(3), 50, player)

class DropofWater(Card):
	def __init__(self, player):
		super().__init__(player)
		self.activated = F_DropofWater

def F_AcidSplash():
	# generate 15 simple acid molecules in the appropriate player's zone
	mp.createMols('SO4', np.zeros(3), 15, player)

class AcidSplash(Card):
	def __init__(self, player):
		super().__init__(player)
		self.activated = F_AcidSplash

def F_DropofLye():
	# generate 15 simple basic molecules in the appropriate player's zone
	mp.createMols('NaOH', np.zeros(3), 15, player)

class DropofLye(Card):
	def __init__(self, player):
		super().__init__(player)
		self.activated = F_DropofLye

def F_ElectronDeathRay():
	# generate a chlorophyll
	mp.createMols('chlorophyll', np.zeros(3), 1, player)

class ElectronDeathRay(Card):
	def __init__(self, player):
		super().__init__()
		self.activated = F_ElectronDeathRay

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

def F_FirinMahLazah():
	uvbeam.on(3)

class FirinMahLazah(Card):
	def __init__(self, player):
		super().__init__()
		self.activated = F_FirinMahLazah

def F_DryingSalts():
	# provide some dehydrating salts
	mp.createMols('MgCO3',0,5,player)
	mp.createMols('CaCO3',0,5,player)
	mp.createMols('MgSO4',0,5,player)

class DryingSalts(Card):
	def __init__(self, player):
		super().__init__()
		self.activated = F_DryingSalts

def F_Thing():
	# do thing
	mp.createMols('thing',0,10,player)

class Thing(Card):
	def __init__(self, player):
		super().__init__()
		self.activated = F_Thing
