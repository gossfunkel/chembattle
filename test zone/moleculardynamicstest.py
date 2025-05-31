from ursina import *
from ursina.shaders import lit_with_shadows_shader
import numpy as np
import math
import sys

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
Entity.rotation_directions  = (1,1,1)

total_T_label				= Text('Total internal heat:',origin=(1.7,-18))
total_T_display				= Text('x',origin=(4,-16))
collisions_count			= 0
#collisions_label			= Text('Collisions:',origin=(-1,-18))
mouseloc_label				= Text('Mouse position:',origin=(-1,-18))
#collisions_display			= Text('x',origin=(-10,-16))
mouseloc_display			= Text('x',origin=(-1,-16))

BC 							= False # periodic boundaries = true; reflective boundaries = false
mass 						= [1.0,16.0]
qs							= [0.41,-0.82] #charges of particles
Temp0 						= 200 # temperature in Kelvin

#min_distance = 0.1
#max_distance = 10

#set_dt = 0.002
#speedlimit = np.array([8,8,8])
#Interaction parameters, roughly copied from
#Xin Bo Zhang et al. Fluid Phase Equilibria 262 (2007) 210–216 doi:10.1016/j.fluid.2007.09.005 
#note that 1kJ/mol = 100 Dal A^2 ps^-2 (for energy conversions from table 1 in that paper)
eps = [[3.24,14.2723],[14.2723,62.87]] # rough parameters for H
sig = [[0.98,2.04845],[2.04845,3.1169]] 
# not sure how these sigma relate to values for van der waals radii i found online - H 1.2, O 1.5 (in Angstroms)
#  could be related to van der waals constant b values in cm^3/mol; H 26.61 O 31.83?
#other constants
Kr  = 148000.0  #spring potential is usually defined as U = (k/2)(r-r_0)^2. I included the /2 here
bl  = 0.9611
Kth = 35300.0  #same explanation as Kr but with bending energy
th0 = 109.47*np.pi/180.0 #bond angle in rad (for water)
#k0 = 1.38064852E-23 	# Boltzmann constant in SI
#k1 = k0* 6.0221409e+26 	# Molar gas constant
#k2 = k1*1.0E20 			# rescale from metres to angstroms
#kb = k2/1.0E24 			# rescale time to picoseconds > Boltzmann constant in A^2 D ps^-2 K^-1
kb  = 0.8314459920816467 # Boltzmann
NA  = 6.0221409e+26 #Avogardos constant x 1000 (g->kg)
ech = 1.60217662E-19 #electron charge in coulombs
kc  = 8.9875517923E9*NA*1E30*ech*ech/1E24 #electrostatic constant in Daltons, electron charges, picosecond, angstrom units

n = 30 # number of atoms
D = 3 # number of spacial dimensions
LL = 15 # max size of system
L = np.zeros([D])+LL

#camera.look_at(atoms[0])
camera.world_position = (LL/2,LL/2,-LL-50)

#def AddAtom():
#	global n
#	n += 1
# 	atomPositions.append(np.random.rand(D))
# 	atomVelocities.append(np.random.rand(D))
# 	atomAccel.append(np.empty(D))
# 	atoms.append(Atom(atomPositions[n, :],atomVelocities[n, :],x,color.random_color()))

#addbutton = Button(origin=(1.5,1.5),text='add atom',scale=0.2,background=color.white,on_click=AddAtom)

class Atom(Entity):
	def __init__(self, pos, vel, x, col=color.green):
		super().__init__(model='sphere',
						world_position=pos,
						scale=0.75,
						color=col,
						collider=None,
						shader=lit_with_shadows_shader)
		self.velocity = vel
		self.index = x

	def update(self):
		self.world_position = atomPositions[self.index, :]
		#self.x 				= self.x - 10
		#self.y 				= self.y - 10
		self.velocity 		= atomVelocities[self.index, :]

# initialise arrays
atomPositions   = (LL*(np.random.rand(n,D))) 		 # initialise all positions randomly
#atomPositions   = LL*(np.empty((n,D))) 		 # initialise all positions randomly
#print(atomPositions)
#atomVelocities  = 10.00*(np.random.rand(n,D)-0.5) # initialise all velocities randomly
atomAccel 		= np.empty((n,D)) # initialise all acceleration to nil
atomVelocities = np.zeros((n,D))				 # initialise all velocities as 0

for x in range(2):
	for y in range(2):
		#print(x)
		e = Entity(model='sphere',color=color.white,scale=2,world_position=((x)*LL,(y)*LL,LL))
		e = Entity(model='sphere',color=color.gray,scale=2,world_position=((x)*LL,(y)*LL,0))
# e = Entity(model='sphere',color=color.white,scale=0.5,world_position=(10,10,-10))
# e = Entity(model='sphere',color=color.white,scale=0.5,world_position=(10,-10,-10))
# e = Entity(model='sphere',color=color.white,scale=0.5,world_position=(-10,-10,-10))
# e = Entity(model='sphere',color=color.white,scale=0.5,world_position=(-10,-10,10))
# e = Entity(model='sphere',color=color.white,scale=0.5,world_position=(-10,10,10))
# e = Entity(model='sphere',color=color.white,scale=0.5,world_position=(-10,10,-10))
# e = Entity(model='sphere',color=color.white,scale=0.5,world_position=(10,-10,10))

#bonds
bnd=[]
for i in range(int(n/3)):
	# structure: [atom,bonded-to,bond-length,gamma(bond-strength)]
	bnd.append([3*i,3*i+1,bl,Kr]) 
	bnd.append([3*i+1,3*i+2,bl,Kr])
bnd=np.array(bnd)

#angles
angs=[]
for i in range(int(n/3)):
	angs.append([3*i,3*i+1,3*i+2,th0,Kth])
angs=np.array(angs)

atoms = []
for x in range(n):
	#if tp[n]:
	#	col = color.red
	#else:
	#	col = color.blue
	atoms.append(Atom(atomPositions[x, :],atomVelocities[x, :],x,color.random_color()))

#Types in groups of three
tp=[0]*n
for i in range(int(n/3)):
    tp[3*i]=0
    atoms[3*i].color = color.red
    tp[3*i+1]=1
    atoms[3*i+2].color = color.green
    tp[3*i+2]=0
    atoms[3*i+1].color = color.blue

#molecule labels
mols=[0]*n
for i in range(int(n/3)):
    mols[3*i]=i
    mols[3*i+1]=i
    mols[3*i+2]=i

#mass and charge arrays
mm=np.array([mass[tp[j]] for j in range(n)])
chrg=np.array([qs[tp[j]] for j in range(n)])

def reflectBC(r,v):
	newv = 1.0 * v
	newr = 1.0 * r
	for i in range(n):
		for j in range(D):
			if newr[i][j]<0:
				newr[i][j]= -newr[i][j]
				newv[i][j]=abs(v[i][j])
			if newr[i][j]>L[j]:
				newr[i][j]= 2.0*L[j]-newr[i][j]
				newv[i][j]=-abs(v[i][j])
	return newr,newv

def LJpot(r,i,sigg,epss):
	sg=np.delete(np.array([sigg[tp[j]] for j in range(n)]),i)
	ep=np.delete(np.array([epss[tp[j]] for j in range(n)]),i)
	for ii in range(n):			#ignore atoms in the same molecule
		if mols[i]==mols[ii]: 
			ep[ii]=0
	drv=r-r[i] #distance in each dimension
	drv=np.delete(drv,i,0) #remove ith element (no self LJ interactions)
	dr=[np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) for a in drv] #absolute distance of that lad
	r6 =(sg/np.array(dr))**6
	r12=(sg/np.array(dr))**12
	LJP=4.0*eps*sum(ep*r6-ep*r12)
	return LJP

#Gradient of Lennard-Jones potential
def dLJp(r,i,sigl,epsl,bdln):
	sg=np.delete(np.array([sigl[tp[j]] for j in range(n)]),i)
	ep=np.array([epsl[tp[j]] for j in range(n)])
	for ii in range(n): #ignore atoms in the same molecule
		if mols[i]==mols[ii]:
			ep[ii]=0
	ep=np.delete(ep,i)
	drv=r-r[i] #distance in each dimension
	drv=np.delete(drv,i,0) #remove ith element (no self LJ interactions)
	dr=[np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) for a in drv] #absolute distance of that lad
	if ((dr[i] < 4) for i in dr):
		print("ep: " + str(ep))
		print("dr: " + str(dr))
		r8 = ep*(sg**6)*(1.0/np.array(dr))**8
		r14=2.0*ep*(sg**12)*(1.0/np.array(dr))**14
		r8v =np.transpose(np.transpose(drv)*r8)
		r14v=np.transpose(np.transpose(drv)*r14)
		r8vs =np.sum(r8v,axis=0)
		r14vs=np.sum(r14v,axis=0)
		dLJP=24.0*(r14vs-r8vs)
	else:
		dLJP=0.0
	return dLJP

#bond length potential
def BEpot(r,bnds):
	bps=np.zeros(n)
	for i in range(n): #loop over all particles
		for j in range(len(bnds)): #check all bonds to see if particle i is bonded
			if bnds[j][0]==i or bnds[j][1]==i:
				if bnds[j][0]==i: #find particle bonded to i
					ii=int(bnds[j][1])
				else:
					ii=int(bnds[j][0])
				dr0=bnds[j][2]
				e0 =bnds[j][3]
				dr=r[i]-r[ii]
				dr2=dr*dr
				adr2=sum(dr2)
				adr=np.sqrt(adr2)
				BE=e0*(adr-dr0)**2
				bps[i]+=BE
	return bps

#gradient of bond length potential (negative force)
def dBEpot(r,bnds):
	bps=np.zeros([n,3])
	for i in range(n): #loop over all particles
		for j in range(len(bnds)): #check all bonds to see if particle i is bonded
			if bnds[j][0]==i or bnds[j][1]==i:
				if bnds[j][0]==i: #find particle bonded to i
					ii=int(bnds[j][1])
				else:
					ii=int(bnds[j][0])
				dr0=bnds[j][2]
				e0 =bnds[j][3]
				dr=r[i]-r[ii]
				dr2=dr*dr
				adr2=sum(dr2)
				adr=np.sqrt(adr2)
				dBE=2.0*e0*(adr-dr0)*dr/adr
				bps[i]+=dBE
	return bps

#gradient of bond angle potential (negative force)
def dBA(r,angs):
	aps=np.zeros([n,3])
	for i in range(n): #loop over all particles
		for j in range(len(angs)): #check all bonds to see if particle i is bonded
			a1=int(angs[j][0])
			a2=int(angs[j][1])
			a3=int(angs[j][2])
			if i==a1 or i==a2 or i==a3:
				th00=angs[j][3] #equilibrium angle
				e0 =angs[j][4] #bending modulus
				if i==a1 or i==a2:
					r1=r[a1]-r[a2] #bond vector 1 (form middle atom to atom 1)
					r2=r[a3]-r[a2] #bond vector 2 (middle atom to atom 2)
				else:
					r1=r[a3]-r[a2] #bond vector 1 (form middle atom to atom 1)
					r2=r[a1]-r[a2] #bond vector 2 (middle atom to atom 2)
				ar1=np.sqrt(sum(r1*r1)) #lengths of bonds
				ar2=np.sqrt(sum(r2*r2))
				dot=sum(r1*r2) #r1 dot r2
				ndot=dot/(ar1*ar2) #normalize dot product by vector lengths i.e. get the cos of angle
				th=math.acos(ndot) #bond angle, theta
				dUdth=-2.0*e0*(th-th00) #-dU/dtheta
				if a1==i or a3==i:
					numerator=(r2/(ar1*ar2))-(dot/(ar1*ar1*ar1*ar2*2.0))
					denominator=np.sqrt(1.0-ndot*ndot)
					dUdr=dUdth*numerator/denominator
					aps[i]+=dUdr
				if i==a2:
					denominator=np.sqrt(1.0-ndot*ndot)
					n1=-(r2+r1)
					n2=dot*r1/(ar1*ar1)
					n3=dot*r2/(ar2*ar2)
					numerator=(n1+n2+n3)/(ar1*ar2)
					dUdr=dUdth*numerator/denominator#
					aps[i]+=dUdr
	return aps

#derivative of coulomb potential (negative force)
def coul(r,i,chrgs):
	q0=chrgs[i]
	qs=1.0*np.array(chrgs)
	for j in range(n):
		if mols[i]==mols[j]:
			qs[j]=0.0
	qs=np.delete(qs,i)
	drv=r-r[i] #distance in each dimension
	drv=np.delete(drv,i,0) #remove ith element (no self LJ interactions)
	dr=[np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) for a in drv] #absolute distance
	r3=q0*qs*kc*((1.0/np.array(dr))**3.0) # Coulomb's law
	FF =np.transpose(np.transpose(drv)*r3) # transpose the distance array, multiply by force, transpose back
	Fs=np.sum(FF,axis=0) # sum the value of axis 0 of the transposed & multiplied array
	#print(Fs)
	return Fs
	#global D
	#return np.zeros(D) disable charge force

def rescaleT(v,T):
	KE = 0.5 * sum(sum(mm*np.transpose(v*v)))
	avKE = KE/n
	Tnow = (2.0/3)*avKE/kb
	lam=np.sqrt(T/Tnow)
	lam=(lam-1.0)*0.5 + 1.0 # uncomment to update temp slowly
	vnew = lam*v
	return vnew

def calculateTemperature(v,n):
	totalvec = np.empty(3)
	for i in range(n):
		totalvec += abs(v[i, :])
	totalvec = np.linalg.norm(totalvec)
	return str(round(totalvec,2))

def updatev(r,v,sigg,epss,a):
	#calculate acceleration:
	F=-np.array([dLJp(r,i,sigg[tp[i]],epss[tp[i]],bnd) for i in range(n)]) #LJ
	F=F-dBEpot(r,bnd) #Bonds
	F=F-dBA(r,angs) #Bond angles
	F=F-np.array([coul(r,i,chrg) for i in range(n)]) #Coulomb
	a=np.transpose(np.transpose(F)/mm) #Force->acceleration
	#update velocity
	#global set_dt
	#newv=v+set_dt*a
	newv=v+time.dt/50*a
	return newv,a

def update():
	global atomPositions
	global atomVelocities
	global atomAccel
	#global time_roller
	#while(time_roller<time.dt):
	atomVelocities, atomAccel = updatev(atomPositions,atomVelocities,sig,eps, atomAccel)
	atomVelocities = rescaleT(atomVelocities,Temp0) #scale to temperature
	atomPositions = atomPositions + atomVelocities * time.dt/50
	# boundaries
	if BC == True:
		atomPositions = atomPositions%L
		atomVelocities = 1.0 * atomVelocities
	else:		
		atomPositions, atomVelocities = reflectBC(atomPositions,atomVelocities)
	#	time_roller += set_dt
	total_T_display.text = calculateTemperature(atomVelocities,n)
	#collisions_display.text = str(collisions_count)
	mouseloc_display.text = str(mouse.collisions)

if __name__ == "__main__":
	app.run()

sys.exit()