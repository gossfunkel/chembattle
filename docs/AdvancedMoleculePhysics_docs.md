AdvancedMoleculePhysics(AMP) docs
===========

### AMP variables
##### n
	(_integer_) number of particles in simulator
##### D
	(_integer_) number of spacial dimensions
##### LL
	(_integer_) scale of the simulation space
##### L
	(_np.ndarray_) LL as a "D"-dimensional array
##### Temp0
	(_integer_) maintained temperature - TODO: CHANGE FOR PRESSURE (walls)
##### atomPositions[]
	(_np.ndarray_) This "D"-dimensional array holds the positions (as vector arrays) of all particles (atoms) in the simulation. It is updated every cycle of the AMP "Update()" loop using the physics procedures.

	Initialised as a Python _list_ before conversion to an _np.ndarray_ - flexible size of list allows populating the array with a dynamic number of molecules. This interconversion to Python _list_s takes place every time the size of this array is modified through the adding or removal of atoms from the simulation. 
##### atomVelocities[]
	(_np.ndarray_) This "D"-dimensional array is identical to the AtomVelocities array in every way, except that it stores the velocity vectors for every particle (atom) in the simulation.
	
	Similarly, it is initialised and converted to a Python _list_ every time its size must be changed.
##### molecules[]
	(_list_) A list of all the _Molecule_ objects currently loaded into the simulation. 
##### tp[]
	(_list_) A list of all the _Atom_ type objects loaded into the simulation. These are our elements/isotopes/ions from which we build our molecules. 
##### a[]
	(_list_) This list will be overwritten every time the AMP "Update()" function runs, with new acceleration vectors for each particle (atom) derived from F=ma - the summed forces on each atom are divided by the mass of the atom to give a new acceleration.
##### Kr
	(_float_) *148000.0* A spring potential for the bond energy, usually defined as U = (k/2)(r-r_0)^2. Can half this value if desired, to taste.
##### Kth
	(_float_) *35300.0*  Same as Kr but for bond bending energy
##### Kb
	(_float_) *0.8314459920816467* The Boltzmann constant
##### NA
	(_float_) *6.0221409e+26* Avogadros constant in *kilograms* (x 1000).
##### ech
	(_float_) *1.60217662E-19* Charge of an electron in *coulombs*.
##### kc
	(_float_) *8.9875517923E9*NA*1E30*ech*ech/1E24* The electrostatic constant in Daltons/electron charges/picosecond/angstrom units.

### CLASS Molecule (inherits from _Entity_)
#### 	__init___(self, nam, position, indx, chldatoms, tpis, velocity=np.zeros(D), bl=[1], th0=109.47*np.pi/180.0)
			Initialises the arrays of all qualities of component atoms, and all collective information.

			The entity visible_self value is initialised to False, making the Molecule itself have no visible form. Instead, the atoms are invisibly packaged by this molecular wrapper, which can be used to handle collisions(/reactions), organise particle data, intelligently look up data with ChemPy, and minimise redundancy. Similarly, the Molecule is not given a model. 

			I may add colliders to the molecules later, to handle the chemical interactions that override the AMP calculations.

			Each molecule has
			n 		= number of atoms in molecule
			index 	= index of molecule in AdvancedMoleculePhysics.molecules[]
			name 	= chempy-compatible chemical name in formulaic form (e.g. 'H20', 'C6H12O6')

			children[] 	= array of all atoms in molecule
			mm[] 		= array of masses of all atoms in molecule
			chrg 		= array of charge values of all atoms in molecule
			covbonds{} 	= dictionary of indices of atoms with (covalent) bonds in molecule to their qualities
							contained array structure: [atom, bonded-to, bond-length, gamma(bond-strength)]
			angles{} 	= dictionary of indices of atoms and the bond angles of their bonds (defaults to float type variable, value '1')
							contained array structure: [prior-atom, atom, following-atom, bond-angle, bending-resistance-force]
			sig{} 		= dictionary of indices of atoms with their sigma constant value
			eps{} 		= dictionary of indices of atoms with their epsilon constant value

#### 	__array__()
			This class overrides "__array__()", allowing it to be treated as an _iterable_. It returns a 2D numpy array. The top level contains three arrays: children (atoms in molecule), mm (masses of children), chrg (charges of children)

#### 	update()
			This procedure is automatically called by Ursina in the game loop. Sets the position and velocity variables of the "Atom" objects in the "children[]" object in this class, taking values from the "atomPositions" and "atomVelocities" np.ndarrays in which calculations are made. This is a likely place for articulation between the render loop and the simulation loop if they are separated for online sync purposes later. 

### AMP procedures:
#### CreateMolecule(location, *mols, player=0, bl=[1], velocity=[np.zeros(D)])
		Initiates a Molecule(name, *location, index, *sig, *eps, *velocity, *atoms), and adds the contents of this molecules to the AMP arrays.

		I will reconstruct this method at some point, because the construction of new molecules is kind of split between this and the actual __init__() function of the Molecule class. Logically, and in terms of cycle efficiency, this seems messy and poorly considered.

		This function is responsible for keeping the numpy arrays properly populated as atoms are added to the simulation. Values are added for Sigma and Epsilon constants per atom, and in order to be utilised for the dLJP function, an average must be calculated between each and the value of every other atom using the 'Lorentz-Bethelot Combining rules'. Each type of atom has its own constants, so for every atom, the of its values in the sigm and epsi arrays
			
#### UpdateMP():
		Called in the game loop. Loops through molecules, and for each, loops through its atoms and calculates the forces on each atom from the "dLJP()", "dBEpot()", "dBA()", and "coul()" procedures. Accelerations of each atom are calculated from these forces, and the "atomVelocities" array is updated using v += a * time.dt/50. "atomPositions" is then updated using r += v * time.dt/50. The game dt value is used for consistent movement across machines, and divided by 50 to scale down movement scale.
