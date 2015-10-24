import interp_sph as interp
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

global gamma,gm1,final_time,dt,alpha,beta
gamma = 1.4
gm1 = gamma-1
final_time = 0.2
dt = 1e-4
nt = int(final_time/dt)
alpha = 1
beta = 1

sound_speed = lambda rho,p: sqrt(gamma*p/rho)

def euler_integrate(u0, rhs):
	return u0 + rhs*dt

def initialise():
	rho_l = 1
	p_l = 1
	u_l = 0

	rho_r = 0.25
	p_r = 0.1
	u_r = 0

	n_particles = 400
	xmin = -0.5
	xmax = 0.5
	dx_l = 0.0015625
	dx_r = 0.00625

	x_l = np.linspace(xmin,0,1+abs(xmin)/dx_l)
	x_r = np.linspace(0,xmax,1+abs(xmin)/dx_r)
	x = np.concatenate((x_l,x_r[1:]))
	
	rho = np.zeros_like(x)
	p = np.zeros_like(x)
	u = np.zeros_like(x)

	rho[x<=0] = rho_l
	rho[x>0] = rho_r 

	p[x<=0] = p_l
	p[x>0] = p_r

	u[:] = 0

	return x, rho, rho*u, p/(gm1*rho)

class Particle:
	def __init__(self, conservative_variables, position):
		self.rho = conservative_variables[0]
		self.rhou = conservative_variables[1]
		self.ener = conservative_variables[2]

		self.u = self.rhou/self.rho
		self.p = (gm1)*self.rho*self.ener
		self.c = sqrt(gamma*self.p/self.rho)
		self.x = position
		self.m = 0.0015625
		self.h = 2*0.00625


def init_particles():
	[x, rho, rhou, ener] = initialise()
	particles = []
	for i in range(len(x)):
		particles.append(Particle([rho[i], rho[i], ener[i]],x[i]))
	return particles

def get_artificial_viscosity(rho_ij, c_ij, r_ij, v_ij):
	if np.dot(vij, rij) < 0:
		mu_ij = h_ij*np.dot(v_ij, r_ij)/(abs(rij)**2 + 1e-6)
		pi_ij = (-alpha*c_ij*mu_ij + beta*mu_ij**2)/rho_ij
		return pi_ij
	else:
		return 0




def euler_rhs(particles):
	n_particles = len(particles)

	rho = np.zeros(n_particles)
	v  = np.zeros_like(rho)
	ener = np.zeros_like(rho)

	for i in range(n_particles):
		for j in range(n_particles):
			r_ij = particles[i].x - particles[j].x
			h_ij = (particles[i].h + particles[j].h)/2
			W_ij = interp.cubic_spline_kernel(r_ij,h_ij)
			W_ij_prime = interp.cubic_spline_derivative(r_ij,h_ij)
			v_ij = (particles[i].v + particles[j].v)/2
			c_i = sound_speed(particles[i].rho, particles[i].p)
			c_j = sound_speed(particles[j].rho, particles[j].p)
			c_ij = np.average([c_j,c_j])
			rho_ij = np.average([particles[i].rho, particles[j].rho])

			pi_ij = get_artificial_viscosity(rho_ij, c_ij, r_ij, v_ij)
			

			rho[i] += particles[i].m*W_ij





