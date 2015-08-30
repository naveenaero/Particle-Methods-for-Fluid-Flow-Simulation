import numpy as np
import cmath
import math
import matplotlib.pyplot as plt

#body parameters
radius = 1
n_body_panels = 4

#boundary parameters
side = 2
n_boundary_panels = side

def distribute_tracers(n_tracers,x_init,y_extreme):
	tracer_initial_location = np.zeros((n_tracers,1),dtype=complex)
	y = np.linspace(-y_extreme,y_extreme,n_tracers)
	for i in range(n_tracers):
		tracer_initial_location[i] = complex(x_init,y[i])
	return tracer_initial_location



def dot_product(v1,v2):
	v2 = v2/np.linalg.norm(v2)
	dot = v1*np.conj(v2)
	return dot.real

class Geometry:
	def __init__(self,n_panels,shape):
		self.n_panels = n_panels
		self.shape = shape

class BodyGeometry(Geometry):

	def __init__(self,radius,n_panels,shape):
		Geometry.__init__(self,n_panels,shape)
		self.radius = radius
		self.nodes = np.zeros((n_panels+1,1),dtype=complex)
		self.cp = np.zeros((n_panels,1),dtype=complex)
		self.normal = np.zeros((n_panels,1),dtype=complex)
		self.cp_theta = np.zeros((n_panels,1))
		
	def set_nodes(self):
		n_nodes = self.n_panels+1
		theta = np.linspace(0,2*np.pi,n_nodes)
		self.nodes = self.radius*np.exp(complex(0,1)*theta)
		self.nodes[-1] = self.nodes[0]
		
	def set_control_points(self):
		n_nodes = self.n_panels+1
		theta = np.linspace(0,2*np.pi,n_nodes)
		for i in range(self.n_panels):
			self.cp_theta[i] = (theta[i]+theta[i+1])/2
			self.cp[i] = (self.nodes[i+1]+self.nodes[i])/2

	def compute_cp_normal(self):
		for i in range(self.n_panels):
			panel_vec = self.nodes[i+1]-self.nodes[i]
			self.normal[i] = self.cp[i]



class BoundaryGeometry(Geometry):
	def __init__(self,side,n_panels,shape):
		Geometry.__init__(self,n_panels,shape)
		self.side = side
		self.shape = shape
		self.n_nodes = (n_panels+1)*4-4
		self.nodes = np.zeros((self.n_nodes,1),dtype=complex)
		self.cp = np.zeros((self.n_nodes,1),dtype=complex)
		self.normal = np.zeros((self.n_nodes,1),dtype=complex)

	def set_nodes(self):
		l = self.side/2.0
		h = self.side/2.0
		panel_length = self.side/self.n_panels
		for i in range(self.n_panels):
			self.nodes[i] = complex(l-i*panel_length , h)
			self.nodes[i+self.n_panels] = complex(-l , h-i*panel_length)
			self.nodes[i+2*self.n_panels] = complex(-l+i*panel_length , -h)
			self.nodes[i+3*self.n_panels] = complex(l , -h+i*panel_length)

	def set_control_points(self):
		for i in range(len(self.nodes)-1):
			self.cp[i] = (self.nodes[i] + self.nodes[i+1])/2
		self.cp[-1] = (self.nodes[0]+self.nodes[-1])/2

	def compute_cp_normal(self):
		for i in range(self.n_panels):
			panel_vec = self.nodes[i+1]-self.nodes[i]
			self.normal[i] = complex(-panel_vec.imag,panel_vec.real)

		
def compute_gamma_coefficients(z,Body,k):
	z1 = Body.nodes[k]
	temp = Body.nodes[k+1]-Body.nodes[k]
	theta = Body.cp_theta[k] + np.pi/2
	l = np.linalg.norm(temp)

	zprime = (z-z1)*cmath.exp(complex(0,-theta))
	
	a1 = (zprime/l)-1
	a2 = a1+1
	b = cmath.log((a1/a2)[0])
	c = cmath.exp(complex(0,theta))
	d = complex(0,-1)/(2*np.pi)
	T1 = np.conj(a1*b + 1)*c*d
	T2 = -np.conj(a2*b + 1)*c*d
	
	return T1,T2

			
			
def compute_geometry_matrix(Body):
	total_panels = Body.n_panels
	A = np.zeros((total_panels+1,total_panels))
	for i in range(total_panels):
		z = Body.cp[i]

		for j in range(total_panels):

			[P1,P2] = compute_gamma_coefficients(z,Body,j)
			gamma_coeff1 = dot_product(P1,Body.cp[i])
			gamma_coeff2 = dot_product(P2,Body.cp[i])
			
			A[i,j] = A[i,j] + gamma_coeff1
			if j+1<=Body.n_panels-1:
				A[i,j+1] = A[i,j+1] + gamma_coeff2
			else:
				A[i,0] = A[i,0] + gamma_coeff2

	A[-1,:] = 1
	return A

Body = BodyGeometry(radius,n_body_panels,'circle')
Body.set_nodes()
Body.set_control_points()
Body.compute_cp_normal()
A = compute_geometry_matrix(Body)
print A
tracer_initial_location = distribute_tracers(10,-2,1.5)


# plt.plot(Body.nodes.real,Body.nodes.imag,linestyle='-')
# plt.plot(Body.cp.real,Body.cp.imag,'o')
# plt.plot(Boundary.nodes.real,Boundary.nodes.imag,linestyle='-')
# plt.plot(Boundary.cp.real,Boundary.cp.imag,'o')

# plt.gca().set_aspect('equal',adjustable='box')
# plt.xlim(-6,6)
# plt.ylim(-6,6)
# plt.show()






