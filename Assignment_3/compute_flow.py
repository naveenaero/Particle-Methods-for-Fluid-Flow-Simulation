import numpy as np
import math
from Geom import *

dt = 0.01
time_horizon = 6
nt = int(time_horizon/dt)
blob_flag = 0

class Solver:
	'''
	class for time integration methods
	'''
	def __init__(self,dt):
		self.dt = dt

	def euler(self,rhs,y0):
		return y0 + rhs*self.dt;

	def rk2(self,updater,DummyTracers,DummyElements,y0,rhs,flag):
		rhs_1 = rhs
		yprime = y0 + rhs_1*self.dt
		if flag==1:
			for i in range(len(rhs)):
				DummyElements[i].CurrentLocation = yprime[i]
			rhs_2 = updater.get_flow_element_velocity(DummyElements)
		else:
			for i in range(len(rhs)):
				DummyTracers[i].CurrentLocation = yprime[i]
			rhs_2 = updater.get_tracer_velocity(DummyTracers,DummyElements) 
		return y0 + (self.dt/2)*(rhs_1 + rhs_2)

class LineVortex:
	def __init__(self,gamma1,gamma2):
		self.gamma1 = gamma1
		self.gamma2 = gamma2

	def get_velocity(self,location,panel_index):
		[coeff_1,coeff_2] = compute_gamma_coefficients(location,Body,panel_index)
		velocity = self.gamma1*coeff_1 + self.gamma2*coeff_2
		return velocity

###########################################################################################
'''
Flow element classes
'''
class Vortex:
	
	def __init__(self,strength,InitialLocation,dx,factor):
		self.strength = strength
		self.InitialLocation = InitialLocation
		self.CurrentLocation = InitialLocation
		self.delta = dx*factor;

	
	def field(self,Location):
 		if blob_flag==1:
 			dist = np.linalg.norm(Location - self.CurrentLocation)
 			dist_vec = Location - self.CurrentLocation 
 			a = complex(-dist_vec.imag,dist_vec.real)
 			b = 1/(2*np.pi*(dist**2+self.delta**2))
 			return np.conj(a*(self.strength*b))
 		else:
 			dist = Location - self.CurrentLocation
 			return complex(0,-1)*self.strength*(1/(2*np.pi*dist))
			
 		    

class Sink:
	def __init__(self,strength,InitialLocation):
		self.strength = strength
		self.InitialLocation = InitialLocation
		self.CurrentLocation = InitialLocation
	def field(self,Location):
		dist = Location - self.CurrentLocation
		return complex(-1,0)*(self.strength/(2*np.pi*dist))

class Source:
	def __init__(self,strength,InitialLocation):
		self.strength = strength
		self.InitialLocation = InitialLocation
		self.CurrentLocation = InitialLocation
	def field(self,Location):
		dist = Location - self.CurrentLocation
		return complex(1,0)*(self.strength/(2*np.pi*dist))

class Doublet:
	def __init__(self,strength,InitialLocation):
		self.strength = strength
		self.InitialLocation = InitialLocation
		self.CurrentLocation = InitialLocation
	def field(self,Location):
		dist = Location - self.CurrentLocation
		return -self.strength/(2*np.pi*dist**2)

class Uniform:
	def __init__(self,strength,direction):
		self.strength = strength
		self.direction = math.radians(direction)
		
	def field(self):
		return self.strength*cmath.e**(-self.direction*complex(0,1))

class Tracer:
	def __init__(self,InitialLocation):
		self.InitialLocation = InitialLocation
		self.CurrentLocation = InitialLocation
		self.velocity = complex(0,0)

###########################################################################################

def compute_total_velocity(FlowElements,TracerElements,Panel):
	velocity_TracerElements = np.zeros((len(TracerElements),1),dtype=complex)
	velocity_FlowElements = np.zeros((len(FlowElements),1),dtype=complex)
	for i in range(len(TracerElements)):
		location = TracerElements[i].CurrentLocation
		for j in range(Body.n_panels):
			temp = Panel[i].get_velocity(location,j)
			velocity_TracerElements[i] = velocity_TracerElements[i] + temp
		velocity_TracerElements[i] = velocity_TracerElements[i] + FlowElements[1].field()

	for i in range(len(FlowElements)):
		location = FlowElements[i].CurrentLocation
		for j in range(Body.n_panels):
			temp = Panel[i].get_velocity(location,j)
			velocity_FlowElements[i] = velocity_FlowElements[i] + temp
		velocity_FlowElements[i] = velocity_FlowElements[i] + FlowElements[1].field()

	return velocity_TracerElements,velocity_FlowElements

###########################################################################################


class UniformFlow:
	def __init__(self,strength,direction):
		self.strength = strength
		self.direction = math.radians(direction)
	def get_velocity(self):
		return self.strength*cmath.exp(-complex(0,1)*self.direction)

def compute_freestream_velocity(location,FlowElements):
	vel = complex(0,0)
	temp = complex(0,0)
	for i in range(Body.n_panels):
		temp = Panel[i].get_velocity(location,i)
		vel = vel + temp
	temp = FlowElements[1].field()
	vel = vel + temp
	return vel



#initiate Flow elements
FlowElements = []
vortex1 = Vortex(1,complex(-1.25,0),0,0)
free_stream = Uniform(1,0)
FlowElements.append(vortex1)
FlowElements.append(free_stream)

#initiate gamma
gamma = np.zeros((Body.n_panels,1))

#initiate solver
sol = Solver(dt)

#initiate tracer location
tracer_current_location = np.zeros((len(tracer_initial_location),1),dtype=complex)

#initiate empty list with tracer type objects
TracerElements = []

#insert tracer objects with initialised tracer locations into empty list
for i in range(len(tracer_initial_location)):
	TracerElements.append(Tracer(tracer_initial_location[i]))

#initialise b in Ax=b
b = np.zeros((Body.n_panels+1,1))

#generate entries in b column vector based on flow elements present
for i in range(Body.n_panels):
	v = FlowElements[1].field()
	b[i] = -1*dot_product(v,Body.cp[i])

#solve for Ax=b system of linear equations using
#least squares estimation to get approximate gamma values
gamma = np.linalg.lstsq(A,b)
gamma = gamma[0]

#initiate a 2D matrix for storing tracer locations at all times
location = np.zeros((nt,len(tracer_initial_location)),dtype=complex)

#initiate an empty list for storing LineVortex elements/objects
Panel = []
#Store Line Vortex elements in list using generated gamma values
for i in range(Body.n_panels):
	if i+1<=Body.n_panels-1:
		Panel.append(LineVortex(gamma[i],gamma[i+1]))
	else:
		Panel.append(LineVortex(gamma[i],gamma[0]))

#Solve for the time evolving system by running it in a loop over time and tracers for position update
for t in range(nt):
	for j in range(len(tracer_initial_location)):
		

		
		TracerElements[j].velocity = compute_freestream_velocity(TracerElements[j].CurrentLocation,FlowElements)
		TracerElements[j].CurrentLocation = sol.euler(TracerElements[j].velocity,TracerElements[j].CurrentLocation)
		tracer_current_location[j] = TracerElements[j].CurrentLocation
		location[t,j] = (tracer_current_location[j])[0]
	

	# plt.plot(Body.nodes.real,Body.nodes.imag,linestyle='-')
	# plt.plot(Body.cp.real,Body.cp.imag,'o')
	# plt.plot(tracer_current_location.real,tracer_current_location.imag,'o')
	# plt.gca().set_aspect('equal',adjustable='box')
	# plt.xlim(-6,6)
	# plt.ylim(-6,6)
	# plt.savefig('./images/'+str(t)+'.png')
	# plt.clf()
	# plt.show()

#plot the object and the path lines
plt.plot(Body.nodes.real,Body.nodes.imag,linestyle='-')
# plt.plot(Body.cp.real,Body.cp.imag,'o')
for j in range(len(tracer_initial_location)):
	plt.plot(location[:,j].real,location[:,j].imag,linestyle='-')
plt.gca().set_aspect('equal',adjustable='box')
plt.xlim(-6,6)
plt.ylim(-6,6)
plt.show()







