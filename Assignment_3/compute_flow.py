import numpy as np
import math
from Geom import *

dt = 0.01
time_horizon = 6
nt = int(time_horizon/dt)

class Solver:
	def __init__(self,dt):
		self.dt = dt
	def euler(self,rhs,y0):
		return y0 + rhs*self.dt

class LineVortex:
	def __init__(self,gamma1,gamma2):
		self.gamma1 = gamma1
		self.gamma2 = gamma2

	def get_velocity(self,location,panel_index):
		[coeff_1,coeff_2] = compute_gamma_coefficients(location,Body,panel_index)
		velocity = self.gamma1*coeff_1 + self.gamma2*coeff_2
		return velocity

class Tracer:
	def __init__(self,InitialLocation):
		self.CurrentLocation = InitialLocation
		self.InitialLocation = InitialLocation
		self.velocity = complex(0,0)

class Freestream:
	def __init__(self,strength,direction):
		self.strength = strength
		self.direction = math.radians(direction)
	def get_velocity(self):
		return self.strength*cmath.exp(-complex(0,1)*self.direction)

def compute_flow_velocity(location,flow):
	vel = complex(0,0)
	temp = complex(0,0)
	for i in range(Body.n_panels):
		temp = Sheet[i].get_velocity(location,i)
		vel = vel + temp
	temp = flow.get_velocity()
	vel = vel + temp
	return vel

#initiate freestream
flow = Freestream(1,0)
gamma = np.zeros((Body.n_panels,1))
sol = Solver(dt)
tracer_currentlocation = np.zeros((len(tracer_location),1),dtype=complex)
tracer = []
for i in range(len(tracer_location)):
	tracer.append(Tracer(tracer_location[i]))

b = np.zeros((Body.n_panels+1,1))
b[-1] = 0


for i in range(Body.n_panels):
	v = flow.get_velocity()
	b[i] = -1*dot_product(v,Body.cp[i])

gamma = np.linalg.lstsq(A,b)
gamma = gamma[0]
# print gamma


location = np.zeros((nt,len(tracer_location)),dtype=complex)

Sheet = []
for i in range(Body.n_panels):
	if i+1<=Body.n_panels-1:
		Sheet.append(LineVortex(gamma[i],gamma[i+1]))
	else:
		Sheet.append(LineVortex(gamma[i],gamma[0]))

for t in range(nt):
	for j in range(len(tracer_location)):
		tracer[j].velocity = compute_flow_velocity(tracer[j].CurrentLocation,flow)
		tracer[j].CurrentLocation = sol.euler(tracer[j].velocity,tracer[j].CurrentLocation)
		tracer_currentlocation[j] = tracer[j].CurrentLocation
		# print tracer_currentlocation[j]
		location[t,j] = (tracer_currentlocation[j])[0]
	

	# plt.plot(Body.nodes.real,Body.nodes.imag,linestyle='-')
	# plt.plot(Body.cp.real,Body.cp.imag,'o')
	# plt.plot(tracer_currentlocation.real,tracer_currentlocation.imag,'o')
	# plt.gca().set_aspect('equal',adjustable='box')
	# plt.xlim(-6,6)
	# plt.ylim(-6,6)
	# plt.savefig('./images/'+str(t)+'.png')
	# plt.clf()
	# plt.show()

plt.plot(Body.nodes.real,Body.nodes.imag,linestyle='-')
plt.plot(Body.cp.real,Body.cp.imag,'o')
for j in range(len(tracer_location)):
	plt.plot(location[:,j].real,location[:,j].imag,linestyle='-')
plt.gca().set_aspect('equal',adjustable='box')
plt.xlim(-6,6)
plt.ylim(-6,6)
plt.show()







