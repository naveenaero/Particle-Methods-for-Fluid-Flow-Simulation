import numpy as np
import cmath
import math
import matplotlib.pyplot as plt

#body parameters
radius = 1
n_body_panels = 50
panel_type = 'linear'

#boundary parameters
side = 2
n_boundary_panels = side



def distribute_tracers(n_tracers, x_init, y_extreme):
	tracer_initial_location = np.zeros((n_tracers),dtype=complex)
	y = np.linspace(-y_extreme,y_extreme,n_tracers)
	for i in range(n_tracers):
		tracer_initial_location[i] = complex(x_init,y[i])
	return tracer_initial_location

def dot_product(v1,v2):
	v2 = v2/np.linalg.norm(v2)
	dot = v1*np.conj(v2)
	return dot.real


class Geometry:
	def __init__(self, n_panels, shape):
		self.n_panels = n_panels
		self.shape = shape
		
class BodyGeometry(Geometry):
	def __init__(self, radius, n_panels, shape, panel_type):
		Geometry.__init__(self,n_panels,shape)
		self.radius = radius
		self.nodes = np.zeros((n_panels+1),dtype=complex)
		self.cp = np.zeros((n_panels),dtype=complex)
		self.normal = np.zeros((n_panels),dtype=complex)
		self.tangent = np.zeros((n_panels),dtype=complex)
		self.cp_theta = np.zeros((n_panels))
		self.panel_type = panel_type
		self.panel_length = complex(0,0)

	def set_nodes(self):
		n_nodes = self.n_panels+1
		theta = np.linspace(0,2*np.pi,n_nodes)
		self.nodes = self.radius*np.exp(complex(0,1)*theta)
		self.nodes[-1] = self.nodes[0]
		self.panel_length = np.linalg.norm(self.nodes[1]-self.nodes[0])
		
	def set_control_points(self):
		n_nodes = self.n_panels+1
		theta = np.linspace(0,2*np.pi,n_nodes)
		for i in range(self.n_panels):
			self.cp_theta[i] = (theta[i]+theta[i+1])/2
			self.cp[i] = (self.nodes[i+1]+self.nodes[i])/2

	def compute_cp_normal(self):
		for i in range(self.n_panels):
			self.normal[i] = self.cp[i]

	def compute_cp_tangent(self):
		for i in range(self.n_panels):
			self.tangent[i] = self.nodes[i+1] - self.nodes[i]




class BoundaryGeometry(Geometry):
	def __init__(self, side, n_panels, shape):
		Geometry.__init__(self,n_panels,shape)
		self.side = side
		self.shape = shape
		self.n_nodes = (n_panels+1)*4-4
		self.nodes = np.zeros((self.n_nodes),dtype=complex)
		self.cp = np.zeros((self.n_nodes),dtype=complex)
		self.normal = np.zeros((self.n_nodes),dtype=complex)

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

		


def compute_linear_gamma_coefficients(z, body, k):
	z1 = body.nodes[k]
	temp = body.nodes[k+1]-body.nodes[k]
	theta = body.cp_theta[k] + np.pi/2
	l = np.linalg.norm(temp)

	zprime = (z-z1)*cmath.exp(complex(0,-theta))
	
	a1 = (zprime/l)-1
	a2 = a1+1
	b = cmath.log((a1/a2))
	c = cmath.exp(complex(0,theta))
	d = complex(0,-1)/(2*np.pi)
	T1 = np.conj(a1*b + 1)*c*d
	T2 = -np.conj(a2*b + 1)*c*d
	
	return T1,T2

def compute_linear_panel_matrix(body):
	total_panels = body.n_panels
	A = np.zeros((total_panels+1,total_panels))

	for i in range(total_panels):
		z = body.cp[i]
		
		for j in range(total_panels):

			[P1,P2] = compute_linear_gamma_coefficients(z,body,j)
			gamma_coeff1 = dot_product(P1,body.cp[i])
			gamma_coeff2 = dot_product(P2,body.cp[i])
			
			A[i,j] = A[i,j] + gamma_coeff1
			if j+1<=body.n_panels-1:
				A[i,j+1] = A[i,j+1] + gamma_coeff2
			else:
				A[i,0] = A[i,0] + gamma_coeff2

	A[-1,:] = 1
	return A



def compute_constant_gamma_coefficient(z, body, k):
	z1 = body.nodes[k]
	temp = body.nodes[k+1]-body.nodes[k]
	theta = body.cp_theta[k] + np.pi/2
	l = np.linalg.norm(temp)

	zprime = (z-z1)*cmath.exp(complex(0,-theta))
	
	a1 = (zprime/l)-1
	a2 = a1+1
	b = cmath.log((a1/a2))
	c = cmath.exp(complex(0,theta))
	d = complex(0,-1)/(2*np.pi)
	T1 = np.conj(b)*c*d
	
	return T1

def compute_constant_panel_matrix(body):
	total_panels = body.n_panels
	A = np.zeros((total_panels+1,total_panels))

	for i in range(total_panels):
		z = body.cp[i]

		for j in range(total_panels):

			P1 = compute_constant_gamma_coefficient(z,body,j)
			gamma_coeff1 = dot_product(P1,body.cp[i])
			
			A[i,j] = A[i,j] + gamma_coeff1
			
	A[-1,:] = 1
	return A



def init_body(panel_type, radius, n_body_panels):
	body = BodyGeometry(radius, n_body_panels, 'circle', panel_type)
	body.set_nodes()
	body.set_control_points()
	body.compute_cp_normal()
	body.compute_cp_tangent()

	if panel_type=='linear':
		A = compute_linear_panel_matrix(body)
	if panel_type=='constant':
		A = compute_constant_panel_matrix(body)
	return body, A





[body, A] = init_body(panel_type, radius, n_body_panels)
# print A








