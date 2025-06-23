import numpy as np 
import scipy as sp 
from scipy.integrate import quad
import matplotlib.pyplot as plt 
import sympy as smp 
import plotly.graph_objects as go 
from IPython.display import HTML

t = smp.symbols('t', positive=True)
x, y, z = smp.symbols('x y z')

r = smp.Matrix([x, y, z])
rPrime = smp.Matrix([smp.cos(4*t), smp.sin(4*t), t])
separationVector = r - rPrime

drPdt = smp.diff(rPrime, t).norm().simplify()
lambd = smp.integrate(drPdt, (t, 0, 2*smp.pi))

integrand = lambd * separationVector / separationVector.norm()**3 * drPdt

dExdt = smp.lambdify([t, x, y, z], integrand[0])
dEydt = smp.lambdify([t, x, y, z], integrand[1])
dEzdt = smp.lambdify([t, x, y, z], integrand[2])

def fieldEnergy(x, y, z):
	return np.array([quad(dExdt, 0, 2*np.pi, args=(x,y,z))[0],
					quad(dEydt, 0, 2*np.pi, args=(x,y,z))[0],
					quad(dEzdt, 0, 2*np.pi, args=(x,y,z))[0]])

x = np.linspace(-2,2,10)
y = np.linspace(-2,2,10)
z = np.linspace(0,2*np.pi,10)
xv, yv, zv = np.meshgrid(x,y,z)

elecField = np.vectorize(fieldEnergy, signature='(),(),()->(n)')(xv,yv,zv)
ex = elecField[:,:,:,0]
ey = elecField[:,:,:,1]
ez = elecField[:,:,:,2]

plt.hist(ex.ravel(), bins=100, histtype='step', label='EM-x')
plt.hist(ey.ravel(), bins=100, histtype='step', label='EM-y')
plt.hist(ez.ravel(), bins=100, histtype='step', label='EM-z')
plt.legend()
plt.xlabel('Electric Field Magnitude')
plt.ylabel('Frequency')
#plt.show()

# round off extreme values
maxValue = 150
ex[ex>maxValue] = maxValue
ey[ey>maxValue] = maxValue
ez[ez>maxValue] = maxValue
ex[ex<-maxValue] = -maxValue
ey[ey<-maxValue] = -maxValue
ez[ez<-maxValue] = -maxValue

# compute line 
linez = np.linspace(0, 2*np.pi, 1000)
linex, liney = np.cos(4*linez), np.sin(4*linez)

# draw plot
data = go.Cone(x=xv.ravel(), y=yv.ravel(), z=zv.ravel(),
				u=ex.ravel(), v=ey.ravel(), w=ez.ravel(),
				colorscale='Inferno', colorbar=dict(title='$x^2$'),
				sizemode='scaled', sizeref=0.5)

layout = go.Layout(title=r'Plot Title',
					scene=dict(xaxis_title=r'x',
						yaxis_title=r'y',
						zaxis_title=r'z',
						aspectratio=dict(x=1,y=1,z=1),
						camera_eye=dict(x=1.2,y=1.2,z=1.2)))

figure = go.Figure(data=data, layout=layout)
figure.add_scatter3d(x=linex,y=liney,z=linez,mode='lines',	line=dict(color='green',width=10))
#figure = plt.figure()

figure.to_html(default_width=1000,default_height=600))
#ax = figure.add_subplot(projection='3d')
#ax.plot_wireframe(ex,ey,ez,rstride=10,cstride=10)
#plt.show()