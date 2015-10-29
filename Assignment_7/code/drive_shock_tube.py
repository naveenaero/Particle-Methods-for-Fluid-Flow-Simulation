import interp_sph as interp
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate
import time

global gamma,gm1,final_time,dt,alpha,beta,start_time
start_time = time.time()
gamma = 1.4
gm1 = gamma-1
final_time = 0.2
dt = 5e-5
nt = int(final_time/dt)
alpha = 1
beta = 1

sound_speed = lambda rho,p: np.sqrt(gamma*p/rho)

def euler_integrate(u0, rhs):
	return u0 + rhs*dt

def initialise():
	'''Initialises the initial condition and the grid
		with ghost particles which are kept dormant
	'''
	rho_l = 1
	p_l = 1
	v_l = 0

	rho_r = 0.25
	p_r = 0.1
	v_r = 0
	

	xmin = -0.5
	xmax = 0.5
	dx_l = 0.0015625
	dx_r = 0.00625
	
	xmin_g = xmin - 0.1
	xmax_g = xmax + 0.1
	
	x_l = np.linspace(xmin_g,0,1+int(abs(xmin_g)/dx_l))
	x_r = np.linspace(0,xmax_g,1+int(abs(xmax_g)/dx_r))
	x = np.concatenate((x_l,x_r[1:]))
	
	ng_l = int(abs(xmin_g - xmin)/dx_l)
	ng_r = int(abs(xmax_g - xmax)/dx_r)

	ng = [ng_l, ng_r]

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
	'''Class for a particle with various particles'''
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
	'''initialises the particles with the initial conditions &
		creates objects for all particles
	'''
	[x, ng, rho, rhou, ener] = initialise()
	particles = []
	for i in range(len(x)):
		particles.append(Particle([rho[i], rhou[i], ener[i]], x[i], ng))
	return particles

def get_properties(particles):
	''' given list of particle objects returns the 
		particle properties
	'''
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

def get_viscosity(rho_ij, c_ij, r_ij, v_ij, h_ij):
	''' computes the numerical viscosity term used
		in the velocity and energy acceleration terms
	'''
	cond1 = (v_ij*r_ij <= 0)
	cond2 = (v_ij*r_ij > 0)

	n_particles = len(rho_ij)
	
	mu_ij = h_ij*(v_ij*r_ij)/(r_ij*r_ij + 1e-5)
	pi_ij = (-alpha*c_ij*mu_ij + beta*mu_ij*mu_ij)/rho_ij

	return cond1*pi_ij

def sound_speed_avg(rho, p, i, j):
	''' computes the sound speed average between
		particle i and j
	'''
	c_i = sound_speed(rho[i], p[i])
	c_j = sound_speed(rho[j], p[j])
	c_ij = np.average([c_i, c_j])
	return c_ij

def euler_rhs(x, rho, v, p, particles ,t):
	''' computes the acceleration terms in the Euler Gas Dynamics
		equations
	'''
	ng_l = particles[0].ng[0]
	ng_r = particles[0].ng[1]

	ng_particles = len(particles)
	n_particles = ng_particles - (ng_l+ng_r)
	
	
	rho_a = np.zeros(n_particles)
	v_a  = np.zeros(n_particles)
	ener_a = np.zeros(n_particles)
	x_sph = np.zeros(n_particles)

	one = np.ones(ng_particles)

	m_i = particles[0].m
	m_j = m_i
	
	r_j = np.copy(x)
	h_ij = particles[0].h*one
	
	v_j = np.copy(v)
	c_j = sound_speed(rho, p)
	
	rho_j = np.copy(rho)
	p_j = np.copy(p)
	
	for i in range(ng_l, ng_l + n_particles):
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
		multiplier_1 = (p_i/(rho_i*rho_i) + p_j/(rho_j*rho_j))

		rho_a[i - ng_l] = np.sum(m_j * W_ij)
		
		# multiple v_a formulations
		v_a[i - ng_l] =  -np.sum(m_j * multiplier * W_ij_prime)
		# v_a[i - ng_l] =  -np.sum(m_j * multiplier_1 * W_ij_prime) - 0.5*np.sum(m_j * pi_ij * W_ij_prime)
		
		# multiple ener_a formulations
		ener_a[i - ng_l] = 0.5 * np.sum(m_j * multiplier * v_ij * W_ij_prime)
		# ener_a[i - ng_l] = np.sum(m_j* (p_i/(rho_i*rho_i)) * v_ij * W_ij_prime) + 0.5*np.sum(m_j * pi_ij * v_ij * W_ij_prime)
		# ener_a[i - ng_l] = np.sum(m_j* multiplier_1* v_ij * W_ij_prime) + 0.5*np.sum(m_j* pi_ij * v_ij * W_ij_prime)
		
		# multiple x_sph formulations
		x_sph[i - ng_l] = 0.5 * np.sum(m_j * (-v_ij) * W_ij /rho_ij)
		# x_sph[i - ng_l] = np.sum(m_j * v_j * W_ij / rho_j)
	
	return x_sph, rho_a, v_a, ener_a


def exactSod(x,t=0.2,x0=0,ql=[1.,1.,0.],qr=[0.25,0.1,0.],gamma=1.4):
	''' Gives exact solution to Sod shock tube problem in order to check for accuracies and compare with computational solution.
	Exact solution is given at x points.
	Algorith taken from http://www.phys.lsu.edu/~tohline/PHYS7412/sod.html
	Ref. Sod, G. A. 1978, Journal of Computational Physics, 27, 1-31.
	'''
	
	#Import stuff
	from sympy.solvers import nsolve
	from sympy import Symbol
	
	#Initiate stuff
	
	shape=x.shape
	x=x.flatten()
	p1 = Symbol('p1')
	[rol,pl,vl]=ql
	[ror,pr,vr]=qr
	
	#Calculate wave velocities and values
	
	cleft=(gamma*pl/rol)**(0.5)
	cright=(gamma*pr/ror)**(0.5)
	m=((gamma-1)/(gamma+1))**0.5
	eq=((2*(gamma**0.5))/(gamma-1)*(1-(p1**((gamma-1)/2/gamma))))-((p1-pr)*(((1-(m**2))**2)*((ror*(p1+(m*m*pr)))**(-1)))**(0.5))
	ppost=float(nsolve(eq,p1,0.))
	rpostByrright=((ppost/pr)+(m*m))/(1+((ppost/pr)*(m*m)))
	vpost=(2*(gamma**0.5))/(gamma-1)*(1-(ppost**((gamma-1)/2/gamma)))
	romid=((ppost/pl)**(1/gamma))*rol
	vshock=vpost*(rpostByrright)/(rpostByrright-1)
	ropost=rpostByrright*ror
	pmid=ppost
	vmid=vpost

	#Calculate locations
	x1=x0-(cleft*t)
	x3=x0+(vpost*t)
	x4=x0+(vshock*t)
	ro=[]
	p=[]
	v=[]
	for i in x:
		csound=((m*m)*(x0-i)/t)+((1-(m*m))*cleft)
		vinst=(1-(m*m))*(((i-x0)/t)+cleft)
		roinst=rol*((csound/cleft)**(2/(gamma-1)))
		pinst=pl*((roinst/rol)**(gamma))
		if i<x1:
			ro.append(rol)
			p.append(pl)
			v.append(vl)
		elif (i>=x4):
			ro.append(ror)
			p.append(pr)
			v.append(vr)
		elif (i<x4) and (i>=x3):
			ro.append(ropost)
			p.append(ppost)
			v.append(vpost)
		elif (i<x3) and (((roinst>rol) and (roinst<romid)) or ((roinst<rol) and (roinst>romid))):
			ro.append(roinst)
			p.append(pinst)
			v.append(vinst)
		else:
			ro.append(romid)
			p.append(pmid)
			v.append(vmid)
			
	#Reshape solutions
	ro=np.array(ro).reshape(shape)
	v=np.array(v).reshape(shape)
	p=np.array(p).reshape(shape)
	
	#calculate conserved variables
	rou = ro*v
	ener=p/(gamma-1)/ro

	return([ro,v,p,rou,ener])


def compute_L2_error(computed, exact):
	''' computes the L2 error norm '''
	return np.linalg.norm(computed-exact)/np.linalg.norm(exact)

def plot_state(x_new, rho_new, v_new, p_new, ener_new):
	''' plots the particle properties '''
	
	exact_sol = exactSod(x_new)
	exact_rho = exact_sol[0]
	exact_v = exact_sol[1]
	exact_p = exact_sol[2]
	exact_ener = exact_sol[4]


	fig = plt.figure(1)

	plt.subplot(2,2,1)
	plt.plot(x_new, rho_new, color='blue',label='computed')
	plt.plot(x_new, exact_rho, color='green',label='analytical')
	plt.title('Density')
	# plt.legend()

	plt.subplot(2,2,2)
	plt.plot(x_new, p_new, color='blue',label='computed')
	plt.plot(x_new, exact_p, color='green',label='analytical')
	plt.title('Pressure')
	# plt.legend()

	plt.subplot(2,2,3)
	plt.plot(x_new, ener_new, color='blue',label='computed')
	plt.plot(x_new, exact_ener, color='green',label='analytical')
	plt.title('Energy')
	# plt.legend()

	plt.subplot(2,2,4)
	plt.plot(x_new, v_new, color='blue',label='computed')
	plt.plot(x_new, exact_v, color='green',label='analytical' )
	plt.title('Velocity')

	# plt.legend()

	e_rho = compute_L2_error(rho_new, exact_rho)
	e_v = compute_L2_error(v_new, exact_v)
	e_ener = compute_L2_error(ener_new, exact_ener)
	e_p = compute_L2_error(p_new, exact_p)
	
	print "Max L2 error = ", max([e_rho, e_v, e_ener, e_p])


	plt.savefig('5.png')
	plt.show()


	

def update_density(x, rho, particles):
	''' function for density summation '''
	ng_l = particles[0].ng[0]
	ng_r = particles[0].ng[1]
	
	ng_particles = len(particles)
	n_particles = ng_particles - (ng_l + ng_r)
	
	one = np.ones(ng_particles)
	m_j = particles[0].m
	r_j = np.copy(x)
	h_ij = one*particles[0].h
	rho_a = np.zeros(n_particles)
	for i in range(ng_l, n_particles + ng_l):
			r_i = x[i]
			r_ij = x[i] - r_j
			W_ij = interp.cubic_spline(r_ij, h_ij)
			rho_a[i - ng_l] = np.sum(m_j*W_ij)
	return rho_a

def integrate(acc, initial_cond, particles):
	''' function for integrating the equations '''
	
	ng_l = particles[0].ng[0]
	ng_r = particles[0].ng[1]
	
	x_sph = acc[0]
	rho_a = acc[1]
	v_a =   acc[2]
	ener_a = acc[3]

	x_new = initial_cond[0]
	rho_new = initial_cond[1]
	v_new = initial_cond[2]
	p_new = initial_cond[3]
	ener_new = initial_cond[4]

	x_new[ng_l:-ng_r] = euler_integrate(x_new[ng_l:-ng_r], x_sph + v_new[ng_l:-ng_r])
	ener_new[ng_l:-ng_r] = euler_integrate(ener_new[ng_l:-ng_r], ener_a)
	v_new[ng_l:-ng_r] = euler_integrate(v_new[ng_l:-ng_r], v_a)
	rho_new[ng_l:-ng_r] = update_density(x_new, rho_new, particles)
	p_new[ng_l:-ng_r] = gm1*ener_new[ng_l:-ng_r]*rho_new[ng_l:-ng_r]
	
	return x_new, rho_new, v_new, p_new, ener_new




def driver(particles, nt, plotter):
	''' function for driving all the other functions '''
	
	[x, rho, v, p] = get_properties(particles)
	
	rho_new = np.copy(rho)
	v_new = np.copy(v)
	p_new = np.copy(p)
	ener_new = p_new/(gm1*rho_new)
	x_new = np.copy(x)
	
	ng_l = particles[0].ng[0]
	ng_r = particles[0].ng[1]
		
	for t in range(nt):
		[x_sph, rho_a, v_a, ener_a] = euler_rhs(x_new, rho_new, v_new, p_new, particles, t)
		
		acc = [x_sph, rho_a, v_a, ener_a]
		initial_cond = [x_new, rho_new, v_new, p_new, ener_new]
		
		[x_new, rho_new, v_new, p_new, ener_new] = integrate(acc, initial_cond, particles)
		
		if (t+1)%nt==0 and plotter:
			plot_state(x_new[ng_l:-ng_r], rho_new[ng_l:-ng_r], v_new[ng_l:-ng_r], p_new[ng_l:-ng_r], ener_new[ng_l:-ng_r])
		print "Time:",(t+1)*dt,"--- %s seconds ---"% (time.time() - start_time)
		
	
def run(nt, plotter):
	''' runs the driver function '''
	particles = init_particles()
	driver(particles, nt, plotter)


run(nt, True)


