import interp_sph as interp
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate
import time

global gamma,gm1,final_time,dt,alpha,beta,start_time
start_time = time.time()
gamma = 1.4
gm1 = gamma-1
final_time = 0.05
dt = 1e-4
nt = int(final_time/dt)
alpha = 1
beta = 1

sound_speed = lambda rho,p: np.sqrt(gamma*p/rho)

def euler_integrate(u0, rhs):
	return u0 + rhs*dt

def initialise():
	rho_l = 1
	p_l = 1
	v_l = 0

	rho_r = 0.25
	p_r = 0.1
	v_r = 0
	
	ng = 0
	n_particles = 400 + 2*ng
	xmin = -0.5
	xmax = 0.5
	dx_l = 0.0015625
	dx_r = 0.00625
	xmin_g = xmin - ng*dx_l
	xmax_g = xmax + ng*dx_r

	x_l = np.linspace(xmin_g,0,1+abs(xmin_g)/dx_l)
	x_r = np.linspace(0,xmax_g,1+abs(xmax_g)/dx_r)
	x = np.concatenate((x_l,x_r[1:]))
	
	rho = np.zeros_like(x)
	p = np.zeros_like(x)
	v = np.zeros_like(x)

	rho[x<=0] = rho_l
	rho[x>0] = rho_r 

	p[x<=0] = p_l
	p[x>0] = p_r

	v[:] = 0

	
	return x, ng, rho, rho*v, p/(gm1*rho)

class Particle:
	def __init__(self, conservative_variables, position, ng):
		self.rho = conservative_variables[0]
		self.rhou = conservative_variables[1]
		self.ener = conservative_variables[2]

		self.v = self.rhou/self.rho
		self.p = gm1*self.rho*self.ener
		self.c = np.sqrt(gamma*self.p/self.rho)
		self.x = position
		self.m = 0.0015625
		self.h = 2*0.00625
		self.ng = ng


def init_particles():
	[x, ng, rho, rhou, ener] = initialise()
	particles = []
	for i in range(len(x)):
		particles.append(Particle([rho[i], rhou[i], ener[i]], x[i], ng))
	return particles

def get_properties(particles):
	n_particles = len(particles)

	x = np.zeros(n_particles)
	rho = np.zeros(n_particles)
	v = np.zeros(n_particles)
	p = np.zeros(n_particles)

	for i in range(n_particles):
		x[i] = particles[i].x
		rho[i] = particles[i].rho
		v[i] = particles[i].v
		p[i] = particles[i].p

	return x, rho, v, p



def get_artificial_viscosity(rho_ij, c_ij, r_ij, v_ij, h_ij):
	if np.dot(v_ij, r_ij) <= 0:
		mu_ij = h_ij*np.dot(v_ij, r_ij)/(abs(r_ij)**2 + 1e-4)
		pi_ij = (-alpha*c_ij*mu_ij + beta*mu_ij**2)/rho_ij
		return pi_ij
	else:
		return 0

def get_viscosity(rho_ij, c_ij, r_ij, v_ij, h_ij):
	cond1 = (v_ij*r_ij <= 0)
	cond2 = (v_ij*r_ij > 0)

	n_particles = len(rho_ij)
	modr = np.abs(r_ij)

	mu_ij = h_ij*(v_ij*r_ij)/(modr*modr + 1e-6)
	pi_ij = (-alpha*c_ij*mu_ij + beta*mu_ij*mu_ij)/rho_ij

	return cond1*pi_ij

def sound_speed_avg(rho, p, i, j):
	c_i = sound_speed(rho[i], p[i])
	c_j = sound_speed(rho[j], p[j])
	c_ij = np.average([c_i, c_j])
	return c_ij

def euler_rhs(x, rho, v, p, particles ,t):
	ng = particles[0].ng
	n_particles = len(x) - 2*ng

	rho_a = np.zeros_like(rho)
	v_a  = np.zeros_like(v)
	ener_a = np.zeros_like(p)
	x_sph = np.zeros_like(x)

	one = np.ones(n_particles)

	m_i = particles[0].m*one
	m_j = np.copy(m_i)
	r_j = np.copy(x)
	h_j = particles[0].h*one
	v_j = np.copy(v)
	c_j = sound_speed(rho, p)
	rho_j = np.copy(rho)
	p_j = np.copy(p)
	h_ij = particles[0].h*one

	# rho = update_density(x, rho, particles)

	for i in range(len(x)):
		r_i  = x[i]
		r_ij = r_i - r_j
		
		W_ij = interp.cubic_spline(r_ij, h_ij)
		W_ij_prime = interp.cubic_spline_kernel_derivative(r_ij, h_ij)
		
		v_i = v[i]*one
		v_ij = v_i - v_j

		c_i = c_j[i]*one
		c_ij = (c_i+c_j)/2
		
		rho_i = rho[i]*one
		rho_ij = (rho_i+rho_j)/2
		
		pi_ij = get_viscosity(rho_ij, c_ij, r_ij, v_ij, h_ij)

		p_i = p[i]*one
		
		multiplier = (p_i/(rho_i*rho_i) + p_j/(rho_j*rho_j) + pi_ij)

		rho_a[i] = np.sum(m_j * W_ij)
		v_a[i] =  -np.sum(m_j * multiplier * W_ij_prime)
		ener_a[i] = 0.5 * np.sum(m_j * multiplier * v_ij * W_ij_prime)
		x_sph[i] = np.sum(m_j * v_j * W_ij / rho_j)
	
	return x_sph, rho_a, v_a, ener_a

def plot_state(x_new, rho_new, v_new, p_new):
	plt.plot(x_new,rho_new)
	plt.show()

def update_density(x, rho, particles):
	n_particles = len(particles)
	rho_a = np.zeros_like(rho)
	for i in range(n_particles):
		for j in range(n_particles):
			r_ij = x[i] - x[j]
			h_ij = (particles[i].h + particles[j].h)/2
			W_ij = interp.cubic_spline_kernel(r_ij, h_ij)
			m_j = particles[j].m
			rho_a[i] += m_j*W_ij
	return rho_a




def driver(particles, nt, plotter):
	
	[x, rho, v, p] = get_properties(particles)
	
	rho_new = np.copy(rho)
	v_new = np.copy(v)
	p_new = np.copy(p)
	ener_new = p_new/(gm1*rho_new)
	x_new = np.copy(x)
	ng = particles[0].ng
	
	for t in range(nt):
		[x_sph, rho_a, v_a, ener_a] = euler_rhs(x_new, rho_new, v_new, p_new, particles, t)
			
		rho_new = rho_a
		v_new = euler_integrate(v_new, v_a)
		ener_new = euler_integrate(ener_new, ener_a)
		x_new = euler_integrate(x_new, x_sph)
		p_new = gm1*ener_new*rho_new
		
		if t%nt==0 and plotter:
			plot_state(x_new, rho_new, v_new, p_new)
		print "Time:",(t+1)*dt,"--- %s seconds ---"% (time.time() - start_time)
		
	
def run(nt, plotter):
	particles = init_particles()
	driver(particles, nt, plotter)

run(nt, False)


