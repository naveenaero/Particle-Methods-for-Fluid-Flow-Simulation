import numpy as np
import math
from Geom import *


dt = 0.05
time_horizon = dt
nt = int(time_horizon/dt)
blob_flag = 0
method = 'doublet_flow'
dr = 0.003 
#rotation rate of cylinder
omega = 0
area = np.pi
circulation = 2*area*omega/(math.sqrt(3))
v_infinity = 1
move_option = 0
doublet_strength = 2*np.pi*radius**2*v_infinity*1
vortex_strength = 0




def error_check_points(Body,dr):
	n_points = Body.n_panels*10
	theta = np.linspace(0,np.pi,n_points)
	check_points = np.zeros((n_points),dtype=complex)
	check_point_velocity = np.zeros((n_points,2))
	check_point_velocity[:,0] = theta
	for i in range(n_points):
		check_points[i] = (Body.radius+dr)*cmath.exp(1j*theta[i])
		
	return check_points,check_point_velocity


[check_points,check_point_velocity] = error_check_points(Body,dr)

def compute_border_velocity(check_points,check_point_velocity,Panel,method):
	if method!='doublet_flow':
		for i in range(len(check_points)):
			vel = complex(0,0)
			for j in range(Body.n_panels):
				vel += Panel[j].get_velocity([check_points[i]],j)
			vel += FlowElements[1].field()
			check_point_velocity[i,1] = np.linalg.norm(vel)
	else:
		for j in range(len(check_points)):
			vel = complex(0,0)
			vel += FlowElements[1].field()
			vel += FlowElements[2].field(check_points[j])
			check_point_velocity[j,1] = np.linalg.norm(vel)
			

	return check_point_velocity
		



'''
Function to convert list elements of position to array
'''
def convert_to_array(FlowElements,TracerElements,option):
	curr_pos = []
	if option=='Elements':
		size = len(FlowElements)
		for i in range(size):
			curr_pos.append(FlowElements[i].CurrentLocation)
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
def convert_from_array(FlowElements,TracerElements,NewLocation,option):
	if option=='Elements':
		size = len(FlowElements)
		for i in range(size):
			FlowElements[i].CurrentLocation = NewLocation[i]
		return FlowElements
	elif option=='Tracer':
		size = len(TracerElements)
		for i in range(size):
			TracerElements[i].CurrentLocation = NewLocation[i]
		return TracerElements



class Solver:
	'''
	class for time integration methods
	'''
	def __init__(self,dt):
		self.dt = dt

	def euler(self,DummyTracers,DummyElements,rhs,y0):
		new_position = []
		new_position.append(y0[0] + rhs[0]*self.dt)
		new_position.append(y0[1]+ rhs[1]*self.dt)

		return new_position[0],new_position[1];

	#flag=0 Tracer and flag=1 FlowELement
	def rk2(self,DummyTracers,DummyFlowElements,rhs,y0,move_option,method):
		#Tracer=0 Element=1
		rhs_1 = []
		rhs_2 = []
		yprime = []
		new_position = []
		
		for j in range(len(rhs)):
			rhs_1.append(rhs[j])
			yprime.append(y0[j]+rhs_1[j]*self.dt)
	
		for i in range(len(DummyTracers)):
			DummyTracers[i].CurrentLocation = (yprime[0])[i]

		if move_option:
			for i in range(len(DummyFlowElements)):
				DummyFlowElements[i].CurrentLocation = (yprime[1])[i]
		
		if method=='Panel':
			solve_linear_system(DummyFlowElements,Panel)
			rhs_2.append(body_local_velocity(DummyTracers,DummyFlowElements,0))
			rhs_2.append(body_local_velocity(DummyTracers,DummyFlowElements,1))
			new_position.append(y0[0] + (self.dt/2)*(rhs_1[0]+rhs_2[0]))
			new_position.append(y0[1] + (self.dt/2)*(rhs_1[1]+rhs_2[1]))
			return new_position[0],new_position[1]
		elif method=="doublet_flow":
			rhs_2.append(local_flow(DummyTracers,DummyFlowElements))
			new_position.append(y0[0] + (self.dt/2)*(rhs_1[0]+rhs_2[0]))
			return new_position[0]
	
		else:
			print "Wrong solution Method used"
		
		
	def rk2_images(self,FlowElements,ImageVortex,ImageVortexNegative,rhs,y0):
		rhs_1 = rhs
		yprime = y0[0] + self.dt*rhs[0]

		FlowElements[0].CurrentLocation = yprime[0]
		
		compute_image_vortex(FlowElements,ImageVortex,ImageVortexNegative)

		v_FlowElement = compute_image_velocity(FlowElements,ImageVortex,ImageVortexNegative)
		
		rhs_2 = v_FlowElement
		out = []
		
		for i in range(len(rhs)):
			out.append(y0[i] + (self.dt/2)*(rhs_1[i] + rhs_2[i]))
		
		return out



class LineVortex:
	def __init__(self,gamma1,gamma2):
		self.gamma1 = gamma1
		self.gamma2 = gamma2

	def get_velocity(self,location,panel_index):
		if Body.panel_type=='linear':
			[coeff_1,coeff_2] = compute_gamma_coefficients(location,Body,panel_index)
			velocity = self.gamma1*coeff_1 + self.gamma2*coeff_2
		elif Body.panel_type=='constant':
			coeff_1 = compute_gamma_coefficients(location,Body,panel_index)
			velocity = self.gamma1*coeff_1
		else:
			print "Error -- Body Panel un-defined"

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
		self.velocity = complex(0,0)

	
	def field(self,Location):
 		if blob_flag==1:
 			dist = np.linalg.norm(Location - self.CurrentLocation)
 			dist_vec = Location - self.CurrentLocation 
 			a = complex(-dist_vec.imag,dist_vec.real)
 			b = 1/(2*np.pi*(dist**2+self.delta**2))
 			return np.conj(a*(self.strength*b))
 		else:
 			dist = Location - self.CurrentLocation
 			return np.conj(complex(0,-1)*self.strength*(1/(2*np.pi*dist)))
			
 		    

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
		return np.conj(-self.strength/(2*np.pi*dist**2))

class Uniform:
	def __init__(self,strength,direction):
		self.strength = strength
		self.direction = math.radians(direction)
		self.CurrentLocation = [complex(0,0)]
		
	def field(self):
		return self.strength*cmath.e**(-self.direction*complex(0,1))

class Tracer:
	def __init__(self,InitialLocation):
		self.InitialLocation = InitialLocation
		self.CurrentLocation = InitialLocation
		self.velocity = complex(0,0)

###########################################################################################


def body_local_velocity(TracerElements,FlowElements,flag):
	temp = complex(0,0)
	if flag==0:
		vel = np.zeros((len(TracerElements),1),dtype=complex)
		for j in range(len(TracerElements)):
			location = TracerElements[j].CurrentLocation
			for i in range(Body.n_panels):
				vel[j] = vel[j] + Panel[i].get_velocity(location,i)
			for i in range(len(FlowElements)):
				if FlowElements[i].__class__.__name__=='Uniform':
					vel[j] = vel[j] + FlowElements[i].field();
				else:
					vel[j] = vel[j] + FlowElements[i].field(location)
		return vel
	else:
		vel = np.zeros((len(FlowElements),1),dtype=complex)
		for j in range(len(FlowElements)):
			if FlowElements[j].__class__.__name__!='Uniform' :
				location = FlowElements[j].CurrentLocation
				for i in range(Body.n_panels):
					vel[j] = vel[j] + Panel[i].get_velocity(location,i)
				for i in range(len(FlowElements)):
					if FlowElements[i].__class__.__name__=='Uniform':
						vel[j] = vel[j] + FlowElements[i].field()
		return vel
		

###########################################################################################
def local_flow(TracerElements,FlowElements):
	velocity = np.zeros((len(TracerElements),1),dtype=complex)
	for i in range(len(velocity)):
		location = TracerElements[i].CurrentLocation
		for j in range(len(FlowElements)):
			if FlowElements[j].__class__.__name__=="Uniform":
				velocity[i] += FlowElements[j].field()
			else:
				velocity[i] += FlowElements[j].field(location)
	return velocity

###########################################################################################

def solve_linear_system(FlowElements,Panel):
	b = np.zeros((Body.n_panels+1,1))
	for i in range(Body.n_panels):
		v = complex(0,0)
		for j in range(len(FlowElements)):
			if FlowElements[j].__class__.__name__=='Uniform':
				v = v + FlowElements[j].field()
			else:
				v = v + FlowElements[j].field(Body.cp[i])

		b[i] = -1*dot_product(v,Body.cp[i])
		b[-1] = circulation
		
	gamma = np.linalg.lstsq(A,b)
	gamma = gamma[0]

	if Body.panel_type=='linear':
		for i in range(Body.n_panels):
			Panel[i].gamma1 = gamma[i]
			if i+1<=Body.n_panels-1:
				Panel[i].gamma2 = gamma[i+1]
			else:
				Panel[i].gamma2 = gamma[0]
	elif Body.panel_type=='constant':
		for i in range(Body.n_panels):
			Panel[i].gamma1 = gamma[i]
			Panel[i].gamma2 = gamma[i]
	else:
		print "Error -- Body Panel un-defined"



###########################################################################################
def compute_image_vortex(FlowElements,ImageVortex,ImageVortexNegative):
	#radial location
	
	image_radial_location = 0
	image_phase = 0
	image_strength = 0
	for i in range(len(FlowElements)):
		if FlowElements[i].__class__.__name__=='Vortex':
			
			location = np.linalg.norm(FlowElements[i].CurrentLocation)
			
			image_strength = -FlowElements[i].strength*radius/location
			image_radial_location = radius**2/location
			
			loc = np.asarray(FlowElements[i].CurrentLocation)
			image_phase = cmath.phase(complex(loc.real,loc.imag))
			 
			image_location = image_radial_location*cmath.exp(1j*image_phase)
			
			ImageVortex.strength = image_strength
			ImageVortex.CurrentLocation = image_location
			
			ImageVortexNegative.strength = -image_strength
			ImageVortexNegative.CurrentLocation = image_location*0

###########################################################################################
def compute_image_velocity(FlowElements,ImageVortex,ImageVortexNegative):
	v_FlowElement = []
	v_ImageVortex = complex(0,0)
	v_ImageVortexNegative = complex(0,0)
	
	for i in range(len(FlowElements)):
		if FlowElements[i].__class__.__name__=='Vortex':
			v_FlowElement.append(ImageVortex.field(FlowElements[i].CurrentLocation))
			v_FlowElement[-1] += ImageVortexNegative.field(FlowElements[i].CurrentLocation)
			
	return v_FlowElement
###########################################################################################


#initiate Flow elements
FlowElements = []
# vortex_location = np.zeros((1),dtype=complex)

vortex1 = Vortex(vortex_strength,[complex(-1.25,0)],0.01,4)
free_stream = Uniform(v_infinity,0)
doublet1 = Doublet(doublet_strength,[complex(0,0)])
FlowElements.append(vortex1)
FlowElements.append(free_stream)
FlowElements.append(doublet1)


#Define Image Vortex and its counterpart at origin
ImageVortex = Vortex(0,complex(0,0),0,0)
ImageVortexNegative = Vortex(0,complex(0,0),0,0)
compute_image_vortex(FlowElements,ImageVortex,ImageVortexNegative)

###########################################################################################

#initiate solver
sol = Solver(dt)

###########################################################################################
#initiate empty list with tracer type objects
TracerElements = []

#insert tracer objects with initialised tracer locations into empty list
for i in range(len(tracer_initial_location)):
	TracerElements.append(Tracer(tracer_initial_location[i]))
###########################################################################################
#initialise b in Ax=b
b = np.zeros((Body.n_panels+1,1))

#generate entries in b column vector based on flow elements present
for i in range(Body.n_panels):
	v = complex(0,0)
	for j in range(len(FlowElements)):
		if FlowElements[j].__class__.__name__=='Uniform':
			v = v + FlowElements[j].field()
		else:
			v = v + FlowElements[j].field(Body.cp[i])
	b[i] = -1*dot_product(v,Body.cp[i])
b[-1] = circulation

#solve for Ax=b system of linear equations using
#least squares estimation to get approximate gamma values
gamma = np.linalg.lstsq(A,b)
gamma = gamma[0]
###########################################################################################

#initiate a 2D matrix for storing tracer anf flow element locations at all times
location_tracer_time = np.zeros((nt,len(tracer_initial_location)),dtype=complex)
location_flowelement_time = np.zeros((nt,1),dtype=complex)
location_vortex_image = np.zeros((nt),dtype=complex)

###########################################################################################

#initiate an empty list for storing LineVortex elements/objects
Panel = []
#Store Line Vortex elements in list using generated gamma values
for i in range(Body.n_panels):
	if Body.panel_type=="linear":
		if i+1<=Body.n_panels-1:
			Panel.append(LineVortex(gamma[i],gamma[i+1]))
		else:
			Panel.append(LineVortex(gamma[i],gamma[0]))
	elif Body.panel_type=="constant":
		Panel.append(LineVortex(gamma[i],gamma[i]))
	else:
		print "Error -- Body Panel un-defined"

###########################################################################################

check_point_velocity = compute_border_velocity(check_points,check_point_velocity,Panel,method)
# print check_point_velocity
###########################################################################################

#Solve for the time evolving system by running it in a loop over time and tracers for position update

if method=='Panel':
	for t in range(nt):
		if t%20==0:
			print(str(int(100.0*t/nt))+"%........")
		
		solve_linear_system(FlowElements,Panel)
		
		# for j in range(len(tracer_initial_location)):
		curr_tracer_position = convert_to_array(FlowElements,TracerElements,'Tracer')
		curr_tracer_velocity = body_local_velocity(TracerElements,FlowElements,0)
		curr_element_position = convert_to_array(FlowElements,TracerElements,'Elements')
		curr_element_velocity = body_local_velocity(TracerElements,FlowElements,1)

		curr_position = [curr_tracer_position,curr_element_position]
		curr_velocity = [curr_tracer_velocity,curr_element_velocity]
		
		[new_tracer_position,new_element_position] = sol.rk2(TracerElements,FlowElements,curr_velocity,curr_position,move_option,method)
		
		TracerElements = convert_from_array(FlowElements,TracerElements,new_tracer_position,'Tracer')
		FlowElements = convert_from_array(FlowElements,TracerElements,new_element_position,'Elements')
		
		for j in range(len(TracerElements)):
			location_tracer_time[t,j] = new_tracer_position[j][0]

		location_flowelement_time[t,0] = new_element_position[0][0]

elif method=='images':
	for t in range(nt):
		if t%20==0:
			print(str(int(100.0*t/nt))+"%........")

		curr_position = [FlowElements[0].CurrentLocation]
		
		v_FlowElement = compute_image_velocity(FlowElements,ImageVortex,ImageVortexNegative)
		curr_velocity = v_FlowElement
		print curr_velocity[0]
		
		new_position = sol.rk2_images(FlowElements,ImageVortex,ImageVortexNegative,curr_velocity,curr_position)

		FlowElements[0].CurrentLocation = new_position[0]
		location_flowelement_time[t,0] = new_position[0][0]
		
		compute_image_vortex(FlowElements,ImageVortex,ImageVortexNegative)

elif method=='doublet_flow':
	for t in range(nt):
		if t%20==0:
			print(str(int(100.0*t/nt))+"%........")

		curr_tracer_position = convert_to_array(FlowElements,TracerElements,'Tracer')
		curr_tracer_velocity = local_flow(TracerElements,FlowElements)
		curr_velocity = [curr_tracer_velocity]
		curr_position = [curr_tracer_position]

		
		new_tracer_position = sol.rk2(TracerElements,FlowElements,curr_velocity,curr_position,move_option,method)


		TracerElements = convert_from_array(FlowElements,TracerElements,new_tracer_position,'Tracer')
		# FlowElements = convert_from_array(FlowElements,TracerElements,new_element_position,'Elements')
		
		for j in range(len(TracerElements)):
			location_tracer_time[t,j] = new_tracer_position[j][0]
		
		



else:
	print "Error -- no method mentioned!"





		


print "[DONE]"


###########################################################################################
def plot_tracers(location_tracer_time,option,file_name):
	if option=='movie':
		theta = np.zeros((4,1))
		radius_vector = np.zeros((4,2),dtype=complex)
		for i in range(len(theta)):
			theta[i] = i*np.pi/2
			
		for t in range(nt):
			if t%5==0:
				plt.plot(Body.nodes.real,Body.nodes.imag,linestyle='-')
				for l in range(len(theta)):
					radius_vector[l,0] = cmath.exp(1j*theta[l][0])
					plt.plot(radius_vector[l,:].real,radius_vector[l,:].imag,linestyle='-')
				
				for j in range(len(location_tracer_time[0,:])):
					if t!=0:
						plt.plot(location_tracer_time[0:t-1,j].real,location_tracer_time[0:t-1,j].imag,linestyle='-')
					else:
						plt.plot(location_tracer_time[0,j].real,location_tracer_time[0,j].imag,linestyle='-')
					plt.plot(location_tracer_time[t,j].real,location_tracer_time[t,j].imag,'o')
				plt.gca().set_aspect('equal',adjustable='box')
				plt.xlim(-6,6)
				plt.ylim(-6,6)
				filename = "./Case_1/%04d.png" % (t/5)
				plt.savefig(filename)
				plt.clf()
			theta = theta + dt*omega

	elif option=='pathlines':
		plt.plot(Body.nodes.real,Body.nodes.imag,linestyle='-')
		for j in range(len(location_tracer_time[0,:])):
			plt.plot(location_tracer_time[:,j].real,location_tracer_time[:,j].imag,linestyle='-')
			plt.plot(location_tracer_time[-1,j].real,location_tracer_time[-1,j].imag,'o')
		plt.gca().set_aspect('equal',adjustable='box')
		plt.xlim(-6,6)
		plt.ylim(-6,6)
		plt.xlabel('X coordinate')
		plt.ylabel('Y coordinate')
		filename = "./"+file_name+".png"
		plt.savefig(filename)
		plt.clf()
	else:
		print "Wrong plotting option"
###########################################################################################

def plot_general(location_flowelement_time,filename):
	plt.plot(Body.nodes.real,Body.nodes.imag,linestyle='-')
	plt.plot(location_flowelement_time[:,0].real,location_flowelement_time[:,0].imag,linestyle='-')

	plt.plot(location_flowelement_time[-1,0].real,location_flowelement_time[-1,0].imag,'o')
	plt.gca().set_aspect('equal',adjustable='box')
	plt.xlim(-6,6)
	plt.ylim(-6,6)
	plt.savefig('./'+filename+'.png')
	plt.show()

def plot_border_velocity(check_point_velocity):
	plt.plot(check_point_velocity[:,0],check_point_velocity[:,1],linestyle='-')
	# plt.gca().set_aspect('equal',adjustable='box')
	# plt.xlim(-2,2)
	# plt.ylim(-1,1)
	plt.show()


plot_border_velocity(check_point_velocity)
# plot_tracers(location_tracer_time,'pathlines','doublet_flow')
# plot_general(location_flowelement_time,'image_method_Vortex')

# location_info = np.zeros((nt,2))
# for i in range(nt):
# 	location_info[i,0] = np.linalg.norm(location_flowelement_time[i,0])
# 	location_info[i,1] = cmath.phase(location_flowelement_time[i,0])
# np.savetxt("images_method_Vortex.csv", location_info, delimiter=",")



np.savetxt("15.csv", check_point_velocity, delimiter=",")

