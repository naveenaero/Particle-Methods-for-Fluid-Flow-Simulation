import numpy as np
import math
from itertools import compress
import time
start_time = time.time()
from Geom import *


dt = 0.1
time_horizon = 1
nt = int(time_horizon/dt)
v_infinity = 1
rho = 1
circulation = 0
mean = 0
Re = 1000.0
viscosity = rho*v_infinity*2*body.radius/Re
sigma = math.sqrt(2*viscosity*dt)




###########################################################################################
############ test functions ################
def test_no_penetration():
    test_velocity = np.zeros(body.n_panels)
    for i in range(body.n_panels):
        test_velocity[i] = dot_product(compute_total_velocity(body, flow_elements, panel, body.cp[i]),body.cp[i])
    print "maximum normal component is: ", max(np.abs(test_velocity))
    
def test_single_vortex_motion(flow_elements):
    flow_elements.append(Vortex(1, complex(-1.25,0), 0.01, 'bo'))
    
###########################################################################################

def get_position_array(flow_elements):
    curr_pos = []
    for i in range(len(flow_elements)):
        curr_pos.append(flow_elements[i].current_location)
    return np.asarray(curr_pos)

def store_positions(flow_elements, current_positions):
    for i in range(len(flow_elements)):
        flow_elements[i].current_location = current_positions[i]

###########################################################################################

class Solver:
    '''
    class for time integration methods
    '''
    def __init__(self, dt):
        self.dt = dt

    def euler(self, flow_elements, panel, rhs, y0):
        return y0 + rhs*self.dt

    def rk2(self, dummy_flow_elements, panel, rhs, y0):
        rhs_1 = rhs
        
        yprime = y0 + rhs_1*self.dt
        
        for i in range(len(flow_elements)):
            dummy_flow_elements[i].current_location = yprime[i]
        
        dummy_freestream_vortex_field = compute_freestream_vortex_field(dummy_flow_elements)

        panel = solve_panel_method(panel, dummy_freestream_vortex_field)
        rhs_2 = compute_flowelement_velocity(dummy_flow_elements, panel)
        
        return y0 + (self.dt/2)*(rhs_1 + rhs_2)
    
###########################################################################################
'''
Flow element classes
'''
class LineVortex:
    def __init__(self, gamma1, gamma2):
        self.gamma1 = gamma1
        self.gamma2 = gamma2

    def get_velocity(self, location, panel_index):
        if body.panel_type=='linear':
            [coeff_1, coeff_2] = compute_linear_gamma_coefficients(location, body, panel_index)
            velocity = self.gamma1*coeff_1 + self.gamma2*coeff_2
        elif body.panel_type=='constant':
            coeff_1 = compute_constant_gamma_coefficient(location, body, panel_index)
            velocity = self.gamma1*coeff_1
        else:
            print "Error -- body panel un-defined"

        return velocity
   
class Vortex:
    
    def __init__(self, strength, initial_location, blob_radius, color):
        self.strength = strength
        self.current_location = initial_location
        self.delta = blob_radius;
        self.velocity = complex(0, 0)
        self.color = color

    
    def blob_field(self, location):
        dist = np.linalg.norm(location - self.current_location)
        dist_vec = (location - self.current_location)
        a = complex(-dist_vec.imag, dist_vec.real)
        b = 1/(2*np.pi*(dist**2+self.delta**2))
        return ((a*(self.strength*b)))

    def chorin_field(self, location):
        dist_vec =  location - self.current_location
        dist = np.linalg.norm(dist_vec)
        a = complex(-dist_vec.imag, dist_vec.real)
        if dist==0:
            return complex(0,0)
        elif dist>self.delta:
            b = 1/(2*np.pi*dist**2)
            return (a*self.strength*b)
        else:
        	b = 1/(2*np.pi*(dist)*self.delta)
        	return (a*self.strength*b)

    def point_vortex_field(self, location):
        dist = location - self.current_location
        return np.conj(complex(0, -1)*self.strength*(1/(2*np.pi*dist)))

class Sink:
    def __init__(self, strength, initial_location):
        self.strength = strength
        self.current_location = initial_location
    def field(self, location):
        dist = location - self.current_location
        return complex(-1, 0)*(self.strength/(2*np.pi*dist))

class Source:
    def __init__(self, strength, initial_location):
        self.strength = strength
        self.current_location = initial_location
    def field(self, location):
        dist = location - self.current_location
        return complex(1, 0)*(self.strength/(2*np.pi*dist))

class Doublet:
    def __init__(self, strength, initial_location):
        self.strength = strength
        self.current_location = initial_location
    def field(self, location):
        dist = location - self.current_location
        return np.conj(-self.strength/(2*np.pi*dist**2))

class Uniform:
    def __init__(self, strength, direction):
        self.strength = strength
        self.direction = math.radians(direction)
        
    def field(self):
        return self.strength*cmath.e**(-self.direction*1j)

class Tracer:
    def __init__(self, initial_location):
        self.current_location = initial_location
        self.velocity = complex(0, 0)

###########################################################################################
def compute_vortex_field(flow_elements, location):
    velocity = complex(0,0)
    for i in range(len(flow_elements)):
        velocity += flow_elements[i].blob_field(location)
    return velocity

def compute_panel_field(panel, location):
    velocity = complex(0,0)
    for i in range(body.n_panels):
        velocity += panel[i].get_velocity(location, i)
    return velocity

###########################################################################################

# def plot_quiver(flow_elements, t):
#     nx,ny = 30,30
#     y,x = np.mgrid[-2:3:nx*1j, -2:3:ny*1j]
#     U = np.zeros((nx,ny));
#     V = np.zeros((nx,ny));
#     for i in range(ny):
#         for j in range(nx):
#             location = complex(x[i,j],y[i,j])
#             temp = compute_total_velocity(body, flow_elements, panel, location)
#             U[i,j] = temp.real
#             V[i,j] = temp.imag


#     plt.plot(body.nodes.real, body.nodes.imag, linestyle='-')
#     plt.quiver(x, y, U, V)
    
#     plt.gca().set_aspect('equal', adjustable='box')
#     plt.xlim(-2, 3)
#     plt.ylim(-2, 3)
#     title = "Velocity field at %f  seconds" % (t*dt)
#     plt.title(title)
#     plt.xlabel(r"X",fontsize=15)
#     plt.ylabel(r"Y",fontsize=15)
#     filename = "./final_case_1/velocity_field/%04d.png" %(t/13)
#     plt.savefig(filename)
#     plt.clf()
    

###########################################################################################
def solve_linear_system(freestream_vortex_field):
    b = np.zeros(body.n_panels+1)
    for i in range(body.n_panels):
        location = (np.abs(body.cp[i])+1e-6)
        v = freestream_vortex_field[i]
        b[i] = -1*dot_product(v, body.cp[i])

    b[-1] = 0
    
    gamma = np.linalg.lstsq(A, b)
    gamma = gamma[0]
    return gamma

def reset_panels(panel, gamma):
    if body.panel_type=='linear':
        for i in range(body.n_panels):
            panel[i].gamma1 = gamma[i]
            if i+1<=body.n_panels-1:
                panel[i].gamma2 = gamma[i+1]
            else:
                panel[i].gamma2 = gamma[0]
    elif body.panel_type=='constant':
        for i in range(body.n_panels):
            panel[i].gamma1 = gamma[i]
            panel[i].gamma2 = gamma[i]
    else:
        print "Error -- body panel un-defined"
    return panel

def init_panels(freestream_vortex_field):
    panel = []
    gamma = solve_linear_system(freestream_vortex_field)
    for i in range(body.n_panels):
        if body.panel_type=="linear":
            if i+1<=body.n_panels-1:
                panel.append(LineVortex(gamma[i], gamma[i+1]))
            else:
                panel.append(LineVortex(gamma[i], gamma[0]))
        elif body.panel_type=="constant":
            panel.append(LineVortex(gamma[i], gamma[i]))
        else:
            print "Error -- body panel un-defined"
    return panel

###########################################################################################

def compute_total_velocity(body, flow_elements, panel, location):
    total_velocity = complex(0,0)
    total_velocity += free_stream.field()
    total_velocity += compute_vortex_field(flow_elements, location)
    total_velocity += compute_panel_field(panel, location)
    return total_velocity

###########################################################################################

def compute_flowelement_velocity(flow_elements, panel):
    flow_element_velocity = np.zeros(len(flow_elements),dtype=complex)
    for i in range(len(flow_elements)):
        location = flow_elements[i].current_location
        flow_element_velocity[i] = compute_total_velocity(body, flow_elements, panel, location) 
    return flow_element_velocity

###########################################################################################

def compute_slip(body, panel, freestream_vortex_field):
    slip_velocity = np.zeros(body.n_panels)
    for i in range(body.n_panels):
        total_velocity = complex(0,0)
        location = (np.abs(body.cp[0])+1e-6)*np.exp(body.cp_theta[i]*1j)
       	total_velocity = freestream_vortex_field[i] + compute_panel_field(panel, location)
        slip_velocity[i] = dot_product(total_velocity,body.tangent[i])
    return slip_velocity

###########################################################################################

def distribute_vortex_blobs(body, panel, flow_elements, slip_velocity):
    blob_radius = body.panel_length/np.pi

    for i in range(body.n_panels):
        new_blobs = []
        blob_quant = 0
        gamma = 0
        remaining_strength = 0
        blob_position = (np.abs(body.cp[i]) + 1e-6 + blob_radius)*np.exp(body.cp_theta[i]*1j)
        
        if slip_velocity[i]<0:
        	gamma_max = -0.1
        	color = 'bo'
        else:
        	gamma_max = 0.1
        	color = 'ro'
        

        blob_quant = int(abs(slip_velocity[i]/gamma_max))
        
        gamma = gamma_max*body.panel_length

        

        for j in range(blob_quant):
            new_blobs.append(Vortex(gamma, blob_position, blob_radius, color))
        # if (slip_velocity[i])-(blob_quant*gamma_max)!=0:
	       #  remaining_strength = (slip_velocity[i] - gamma_max*blob_quant)*body.panel_length
	       #  new_blobs.append(Vortex(remaining_strength, blob_position, blob_radius, color))
        
        flow_elements += new_blobs

###########################################################################################

def diffuse(flow_elements):
    current_positions = get_position_array(flow_elements)
    current_positions.real += np.random.randn(len(current_positions))*sigma + mean
    current_positions.imag += np.random.randn(len(current_positions))*sigma + mean
    current_positions = mirror(current_positions)
    store_positions(flow_elements, current_positions)

###########################################################################################

def mirror(current_positions):
	reflection_radius = np.abs(body.radius+1e-6)
	check_location = np.abs(current_positions)<reflection_radius
	particles_outside = list(compress(xrange(len(check_location)), check_location))
	if len(particles_outside)!=0:
		for item in particles_outside:
			old_radius = np.abs(current_positions[item])
			old_theta = cmath.phase(current_positions[item])
			new_radius = old_radius + 2*(reflection_radius-old_radius)
			new_theta = old_theta
			current_positions[item] = new_radius*cmath.exp(new_theta*1j)
	return current_positions
   
###########################################################################################

def compute_freestream_vortex_field(flow_elements):
	other_velocity = np.zeros(body.n_panels,dtype=complex)
	for i in range(body.n_panels):
		location = (np.abs(body.cp[i])+1e-6)*np.exp(body.cp_theta[i]*1j)
		other_velocity[i] = compute_vortex_field(flow_elements, location)
		other_velocity[i] += free_stream.field()
	return other_velocity

###########################################################################################

# def plot_particles(flow_elements, t, file_factor):
#     plt.plot(body.nodes.real, body.nodes.imag,'k-')
#     for i in range(len(flow_elements)):
#         plt.plot(flow_elements[i].current_location.real, flow_elements[i].current_location.imag,flow_elements[i].color,markersize=2)

#     plt.gca().set_aspect('equal', adjustable='box')
#     plt.xlim(-2, 6)
#     plt.ylim(-2, 2)
#     title = "Particle positions at %f seconds" % (t*dt)
#     plt.title(title)
#     plt.xlabel(r"X",fontsize=15)
#     plt.ylabel(r"Y",fontsize=15)
#     filename = "./final_case_1/particles/%04d.png" % (t/file_factor)
#     plt.savefig(filename)
#     plt.clf()

###########################################################################################

# def plot_drag(drag):
# 	time = np.zeros(len(drag))
# 	D = np.zeros(len(drag))
# 	for i in range(len(drag)):
# 		time[i] = drag[i][0]
# 		D[i] = drag[i][1]
# 	plt.plot(time,D)
# 	plt.xlabel(r"$Time(s)$",fontsize=12)
# 	plt.ylabel(r"$C_d$",fontsize=15)
# 	plt.title(r"$Drag\ Coefficient(C_d)\ vs\ Time$")
# 	plt.savefig("./final_case_1/Drag.png")

###########################################################################################


def solve_panel_method(panel, freestream_vortex_field):
    gamma = solve_linear_system(freestream_vortex_field)
    panel = reset_panels(panel, gamma)
    return panel

###########################################################################################

def convect(flow_elements, panel):
    flow_element_velocity = compute_flowelement_velocity(flow_elements, panel)
    flow_element_position = get_position_array(flow_elements)
    new_flow_element_position = sol.rk2(flow_elements, panel, flow_element_velocity, flow_element_position)
    store_positions(flow_elements, new_flow_element_position)

###########################################################################################

def compute_vortex_momentum(flow_elements, vortex_momentum, t):
	particle_momentum = complex(0,0)
	for vortex in flow_elements:
		particle_momentum += vortex.strength*complex(vortex.current_location.imag,-vortex.current_location.real)
	particle_momentum *= rho
	vortex_momentum.append([t*dt,particle_momentum])

###########################################################################################

def compute_drag_force(vortex_momentum):
	new_momentum = []
	n_avg = 1
	for i in range(0,len(vortex_momentum),n_avg):
		momentum = complex(0,0)
		time = 0
		for j in range(i,i+n_avg):
			momentum += vortex_momentum[j][1]
			time += vortex_momentum[j][0]
		momentum /= n_avg
		time /= n_avg
		new_momentum.append([time,momentum])
	drag = []
	for i in range(len(new_momentum)-1):
		dt = new_momentum[i+1][0] - new_momentum[i][0]
		di = (new_momentum[i+1][1] - new_momentum[i][1])
		drag.append([new_momentum[i][0],abs(-di/dt)])
	return drag
	

###########################################################################################
def solve_flow(flow_elements, panel):
    for t in range(nt):
        if t%10==0:
            print(str(int(100.0*t/nt))+"%........")


        print "###############################"
        print "Time:",t*dt,", Number of Blobs:",len(flow_elements)
        print "###############################"
        
        freestream_vortex_field = compute_freestream_vortex_field(flow_elements)
        
        #Solve Panel Method for no-penetration
        panel = solve_panel_method(panel, freestream_vortex_field)

        
        #compute and store slip velocity at the control points
        slip_velocity = compute_slip(body, panel, freestream_vortex_field)
        
        
        #plot the quiver plot of velocity
        # if t==nt-1:
        # 	plot_quiver(flow_elements, t)

        
        #convect the flow elements
        convect(flow_elements, panel)
        
		#distribute vortex blobs at the control points
        distribute_vortex_blobs(body, panel, flow_elements, slip_velocity)
        
        
        #Diffuse all the flow elements
        diffuse(flow_elements)

        
        # plot the particles at every time step
        # file_factor = 1
       	# if (t+1)%file_factor==0:
	       #  plot_particles(flow_elements, t, file_factor)

		#Calculate vortex momentum
       	compute_vortex_momentum(flow_elements, vortex_momentum, t)
       	
    print "[DONE]"
    
###########################################################################################

#initiate Flow elements
flow_elements = []
vortex_momentum = []
free_stream = Uniform(v_infinity,0)


###########################################################################################

#initiate solver
sol = Solver(dt)

###########################################################################################

#initiate panels
freestream_vortex_field = compute_freestream_vortex_field(flow_elements)
panel = init_panels(freestream_vortex_field)

###########################################################################################
#Solve for the time evolving system
solve_flow(flow_elements, panel)

###########################################################################################
#Compute and plot drag
drag = compute_drag_force(vortex_momentum)
# plot_drag(drag)
print("--- %s seconds ---" % (time.time() - start_time))

###########################################################################################










