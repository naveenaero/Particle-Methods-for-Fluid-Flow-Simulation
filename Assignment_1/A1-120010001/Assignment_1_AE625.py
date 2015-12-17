import cmath as cmath
import numpy as np
import matplotlib.pyplot as plt
import math as math


option  = input('Want the Flow Element to move(1) or not move(0)?:')
method = ['euler']

time_horizon = 4
dt0 = 0.01
dt = np.array([[dt0]])#[dt0/2],[dt0/4],[dt0/8]])
	
error = np.zeros((len(dt),2))
num = 0



'''
Function to compute analytical solution
'''
def ComputeAnalyticalSolution(time_horizon,dt,Elements):
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
def ConvertToArray(Elements,TracerElements,option):
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
def ConvertFromArray(Elements,TracerElements,NewLocation,option):
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
def StorePositionData(Elements,ElementPosition,index):
	size = len(Elements)
	for i in range(size):
		if Elements[i].kind!='uniform':
			ElementPosition[index,i] = Elements[i].CurrentLocation
	return ElementPosition

'''
Function for storing position data of Tracer Elements
'''
def StoreTracerData(TracerElements,TracerPosition,index):
	size = len(TracerElements)
	for i in range(size):
		TracerPosition[index,i] = TracerElements[i].CurrentLocation
	return TracerPosition


class solver:
	'''
	class for time integration methods
	'''
	def __init__(self,dt):
		self.dt = dt

	def Euler(self,y0,RHS):
		return y0 + RHS*self.dt;

	def RK2(self,DummyTracers,DummyElements,y0,RHS,flag):
		RHS_1 = RHS
		yprime = y0 + RHS_1*self.dt
		if flag==1:
			for i in range(len(RHS)):
				DummyElements[i].CurrentLocation = yprime[i]
			RHS_2 = updater.GetFlowElementVelocity(DummyElements)
		else:
			for i in range(len(RHS)):
				DummyTracers[i].CurrentLocation = yprime[i]
			RHS_2 = updater.GetTracerVelocity(DummyTracers,DummyElements) 
		return y0 + (self.dt/2)*(RHS_1 + RHS_2)

class FlowUpdate:
	'''
	class for updating the flow and getting velocity of elements
	'''
	def __init__(self,option):
		self.option = option;
	
	def UpdatePosition(self,DummyTracers,DummyElements,CurrentPos,RHS,method,flag):
		if method=='euler':
			CurrentPos = sol.Euler(CurrentPos,RHS);
		elif method=='RK2':
			CurrentPos = sol.RK2(DummyTracers,DummyElements,CurrentPos,RHS,flag);
		else:
			print('Please enter correct solver type!!')
		return CurrentPos

	def GetFlowElementVelocity(self,Elements):
		size = len(Elements)
		velocity = np.zeros(size)*complex(0,0)
		if self.option==1:
			for j in range(size):
				VelocityIncrement = complex(0,0)
				Location = Elements[j].CurrentLocation
				if Elements[j].kind!='uniform':
					for i in range(size):
						if i!=j:
							VelocityIncrement = Elements[i].Field(Location)
							velocity[j] = velocity[j] + VelocityIncrement
			return np.conj(velocity)
				
		else:
			return np.conj(velocity)
		

	def GetTracerVelocity(self,TracerElements,Elements):
		size_TracerElem = len(TracerElements)
		size_Elements = len(Elements)
		velocity = np.zeros(size_TracerElem)*complex(0,0)

		for j in range(size_TracerElem):
			Location = TracerElements[j].CurrentLocation
			VelocityIncrement = complex(0,0)
			for i in range(size_Elements):
				VelocityIncrement = Elements[i].Field(Location)
				velocity[j] = velocity[j] + VelocityIncrement
		return np.conj(velocity)

'''
Flow element classes
'''
class Vortex:
	kind = 'vortex'
	def __init__(self,strength,InitialLocation):
		self.strength = strength
		self.InitialLocation = InitialLocation
		self.CurrentLocation = InitialLocation

	
	def Field(self,Location):
		dist = Location - self.CurrentLocation
		return complex(0,-1)*(self.strength/(2*np.pi*dist))

class Sink:
	kind = 'sink'
	def __init__(self,strength,InitialLocation):
		self.strength = strength
		self.InitialLocation = InitialLocation
		self.CurrentLocation = InitialLocation
	def Field(self,Location):
		dist = Location - self.CurrentLocation
		return complex(-1,0)*(self.strength/(2*np.pi*dist))

class Source:
	kind = 'source'
	def __init__(self,strength,InitialLocation):
		self.strength = strength
		self.InitialLocation = InitialLocation
		self.CurrentLocation = InitialLocation
	def Field(self,Location):
		dist = Location - self.CurrentLocation
		return complex(1,0)*(self.strength/(2*np.pi*dist))

class Doublet:
	kind = 'doublet'
	def __init__(self,strength,InitialLocation):
		self.strength = strength
		self.InitialLocation = InitialLocation
		self.CurrentLocation = InitialLocation
	def Field(self,Location):
		dist = Location - self.CurrentLocation
		return -self.strength/(2*np.pi*dist**2)

class Uniform:
	kind = 'uniform'
	def __init__(self,strength,direction):
		self.strength = strength
		self.direction = math.radians(direction)
		self.CurrentLocation = 0
	def Field(self,CurrentLocation):
		return self.strength*cmath.e**(-self.direction*complex(0,1))

class Tracer:
	kind = 'tracer'
	def __init__(self,InitialLocation):
		self.InitialLocation = InitialLocation
		self.CurrentLocation = InitialLocation

'''
Function to read flow element data from txt file
'''
def GetFlowElements():
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
				Elements.append(Vortex(strength,pos))
			elif line[0]=='sink':
				Elements.append(Sink(strength,pos))
			elif line[0]=='source':
				Elements.append(Source(strength,pos))
			elif line[0]=='doublet':
				Elements.append(Doublet(strength,pos))
			else:
				prin("Error! no flow element found")

	return TracerElements,Elements



'''
code snippet to get the solution for various schemes and for different time step size
'''
for item in method:
	solution_method = item 

	for u in range(len(dt)):
		[TracerElements,Elements] = GetFlowElements()
		NumberOfTracers = len(TracerElements)
		NumberOfElements = len(Elements)

		# AnalayticalSolution = ComputeAnalyticalSolution(time_horizon,dt[u],Elements)

		nt = int(time_horizon/dt[u])
		print nt,dt[u]

		sol = solver(dt[u])
		updater = FlowUpdate(option)
		t = np.linspace(0,time_horizon,nt+1)


		ElementPosition = np.zeros((nt+1,NumberOfElements))*complex(0,0)
		ElementPosition = StorePositionData(Elements,ElementPosition,0)


		TracerPosition = np.zeros((nt+1,NumberOfTracers))*complex(0,0)
		TracerPosition = StoreTracerData(TracerElements,TracerPosition,0)

		'''
		Loop for Time Integration
		'''
		for k in range(nt):
			# Solve for Tracer motion
			if NumberOfTracers>0:
				TracerVelocity = updater.GetTracerVelocity(TracerElements,Elements)
				curr_pos = ConvertToArray(Elements,TracerElements,'Tracer')
				NewLocation = updater.UpdatePosition(TracerElements,Elements,curr_pos,TracerVelocity,solution_method,0)
				TracerElements = ConvertFromArray(Elements,TracerElements,NewLocation,'Tracer')
				TracerPosition = StoreTracerData(TracerElements,TracerPosition,k+1)
			# SOlve for Flow ELement motion
			ElemVelocity = updater.GetFlowElementVelocity(Elements)
			curr_pos = ConvertToArray(Elements,TracerElements,'Elements')
			NewLocation = updater.UpdatePosition(TracerElements,Elements,curr_pos,ElemVelocity,solution_method,1)
			Elements = ConvertFromArray(Elements,TracerElements,NewLocation,'Elements')
			ElementPosition = StorePositionData(Elements,ElementPosition,k+1)

			
		# calculate the relative error
		# error[u,num] = np.linalg.norm(abs(AnalayticalSolution[-1]-ElementPosition[-1,1]))/np.linalg.norm(abs(AnalayticalSolution[-1]))

	num = num+1


print error
'''
code snippet for plotting the error convergence plots
'''
################################################################################
# fig, ax = plt.subplots()
# ax.loglog(dt,error[:,0],label='euler')
# ax.loglog(dt,error[:,1],label='RK2')
# legend = ax.legend(loc='upper left')

# # The frame is matplotlib.patches.Rectangle instance surrounding the legend.
# frame = legend.get_frame()
# frame.set_facecolor('0.90')

# # Set the fontsize
# for label in legend.get_texts():
#     label.set_fontsize('large')

# for label in legend.get_lines():
#     label.set_linewidth(1.5)  # the legend line width
# plt.xlabel('Log(dt)')
# plt.ylabel('Log(error)')
# slope_euler,intercept_euler = np.polyfit(dt,error[:,0],1)
# slope_RK2,intercept_RK2 = np.polyfit(dt,error[:,1],1)
# print ("Euler_slope="),slope_euler
# print ("RK2_slope="),slope_RK2
# plt.savefig('convergence.png')
###############################################################################



##############################################################################
# '''
# code snippet for plotting Flow element position history
# '''
# for j in range(NumberOfElements):
# 	plt.plot(ElementPosition[:,j].real,ElementPosition[:,j].imag)
# 	plt.plot(ElementPosition[-1,j].real,ElementPosition[-1,j].imag,'ro')
# plt.xlim(-4, 4)
# plt.ylim(-2, 2)
# plt.xlabel('X coordinate')
# plt.ylabel('Y coordinate')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.title('Flow Elements')
# # plt.savefig('vortex_couple_euler.png')
# plt.show()
# plt.clf()

'''
code snippet for plotting tracer element position history
'''
for j in range(NumberOfTracers):
	plt.plot(TracerPosition[:,j].real,TracerPosition[:,j].imag)
	plt.plot(TracerPosition[-1,j].real,TracerPosition[-1,j].imag,'ro')
for j in range(NumberOfElements):
	if Elements[j].kind!='uniform':
		plt.plot(ElementPosition[0,j].real,ElementPosition[0,j].imag,'ro')
plt.xlim(-4, 4)
plt.ylim(-2, 2)
plt.xlabel('X coordinate')
plt.ylabel('Y coordinate')
plt.gca().set_aspect('equal', adjustable='box')
plt.title('Tracer + Flow Elements')
#plt.savefig('source_sink_tracer_new.png')
plt.show()
plt.clf()






