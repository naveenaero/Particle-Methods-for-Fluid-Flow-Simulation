import interp_sph as interp
import numpy as np
import matplotlib.pyplot as plt

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
	np_l = 320
	np_r = n_particles - np_l
	dx_l = 0.0015625
	dx_r = 0.00625

	x_l = np.linspace(xmin,0,1+abs(xmin)/dx_l)
	x_r = np.linspace(0,xmax,1+abs(xmin)/dx_r)
	x = np.concatenate((x_l,x_r[1:]))
	


initialise()

