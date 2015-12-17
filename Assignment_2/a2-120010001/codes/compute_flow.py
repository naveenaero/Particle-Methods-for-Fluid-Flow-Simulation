import cmath as cmath
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import interpolate
from write_flow_elements import dx

# delta/dx factor; set around 3 for good results
delta_dx_scaling_factor = 5
# blob_flag=1 ==> use vortex blob
# blob_flag=0 ==> do not use vortex blob
blob_flag = 1


'''
Function to compute analytical solution
'''
def compute_analytical_solution(time_horizon,dt,Elements):
	radius = Elements[0].InitialLocation - Elements[1].InitialLocation
	radius = math.sqrt(radius.real**2 + radius.imag**2)/2
	v = Elements[0].strength/(2*np.pi*2*radius)
	omega = v/radius

	t = np.linspace(0,time_horizon,(time_horizon/dt)+1)
	position = np.zeros((len(t),1))*complex(0,0)
	for i in range(len(t)):
		position[i] = complex(radius*math.cos(omega*t[i]),radius*math.sin(omega*t[i]))
	return position

'''
Function to convert list elements of position to array
'''
def convert_to_array(Elements,TracerElements,option):
	curr_pos = []
	if option=='Elements':
		size = len(Elements)
		for i in range(size):
			curr_pos.append(Elements[i].CurrentLocation)
		return np.asarray(curr_pos)
	elif option=='Tracer':
		size = len(TracerElements)
		for i in range(size):
			curr_pos.append(TracerElements[i].CurrentLocation)
		return np.asarray(curr_pos)
	else:
		print("Error, Incorrect input option")


'''
Function to convert array position elements to list
'''
def convert_from_array(Elements,TracerElements,NewLocation,option):
	if option=='Elements':
		size = len(Elements)
		for i in range(size):
			Elements[i].CurrentLocation = NewLocation[i]
		return Elements
	elif option=='Tracer':
		size = len(TracerElements)
		for i in range(size):
			TracerElements[i].CurrentLocation = NewLocation[i]
		return TracerElements

'''
Function for storing position data of Flow Elements
'''
def store_position_data(Elements,ElementPosition,index):
	size = len(Elements)
	for i in range(size):
		if Elements[i].__class__.__name__!='Uniform':
			ElementPosition[index,i] = Elements[i].CurrentLocation
	return ElementPosition

'''
Function for storing position data of Tracer Elements
'''
def store_tracer_data(TracerElements,TracerPosition,index):
	size = len(TracerElements)
	for i in range(size):
		TracerPosition[index,i] = TracerElements[i].CurrentLocation
	return TracerPosition


class Solver:
	'''
	class for time integration methods
	'''
	def __init__(self,dt):
		self.dt = dt

	def euler(self,y0,rhs):
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

class FlowUpdate:
	'''
	class for updating the flow and getting velocity of elements
	'''
	def __init__(self,option):
		self.option = option;
	
	def update_position(self,sol,updater,DummyTracers,DummyElements,CurrentPos,rhs,method,flag):
		if method=='euler':
			CurrentPos = sol.euler(CurrentPos,rhs);
		elif method=='rk2':
			CurrentPos = sol.rk2(updater,DummyTracers,DummyElements,CurrentPos,rhs,flag);
		else:
			print('Please enter correct Solver type!!')
		return CurrentPos

	def get_flow_element_velocity(self,Elements):
		size = len(Elements)
		velocity = np.zeros(size)*complex(0,0)
		if self.option==1:
			for j in range(size):
				VelocityIncrement = complex(0,0)
				Location = Elements[j].CurrentLocation
				if Elements[j].__class__.__name__!='Uniform':
					for i in range(size):
						if i!=j:
							VelocityIncrement = Elements[i].field(Location)
							velocity[j] = velocity[j] + VelocityIncrement
			return np.conj(velocity)
				
		else:
			return np.conj(velocity)
		

	def get_tracer_velocity(self,TracerElements,Elements):
		size_TracerElem = len(TracerElements)
		size_Elements = len(Elements)
		velocity = np.zeros(size_TracerElem)*complex(0,0)

		for j in range(size_TracerElem):
			Location = TracerElements[j].CurrentLocation
			VelocityIncrement = complex(0,0)
			for i in range(size_Elements):
				VelocityIncrement = Elements[i].field(Location)
				velocity[j] = velocity[j] + VelocityIncrement
		return np.conj(velocity)

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
		self.CurrentLocation = 0
	def field(self,CurrentLocation):
		return self.strength*cmath.e**(-self.direction*complex(0,1))

class Tracer:
	def __init__(self,InitialLocation):
		self.InitialLocation = InitialLocation
		self.CurrentLocation = InitialLocation

'''
Function to read flow element data from txt file
'''
def get_flow_elements():
	f = open('flow_elements.txt','r')
	NumberOfLines = sum(1 for _ in f)
	f.seek(0,0)
	TracerElements = []
	Elements = []
	for i in range(NumberOfLines):
		line = (f.readline()).split()
		if line[0]=='tracer':
			x = float(line[1])
			y = float(line[2])
			pos = complex(x,y)
			tracer_elem = Tracer(pos)
			TracerElements.append(tracer_elem)
		
		elif line[0]=='uniform':
			strength = float(line[1])
			direction = float(line[2])
			Elements.append(Uniform(strength,direction))

		else:
			x = float(line[1])
			y = float(line[2])
			pos = complex(x,y)
			strength = float(line[3])
			if line[0]=='vortex':
				Elements.append(Vortex(strength,pos,dx,delta_dx_scaling_factor))
			elif line[0]=='sink':
				Elements.append(Sink(strength,pos))
			elif line[0]=='source':
				Elements.append(Source(strength,pos))
			elif line[0]=='doublet':
				Elements.append(Doublet(strength,pos))
			else:
				prin("Error! no flow element found")

	return TracerElements,Elements









