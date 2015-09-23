import numpy as np
import math
from itertools import compress
from Geom import *


dt = 0.075
time_horizon = 3
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
    # flow_elements.append(Vortex(1, complex(-1.5,0), 0.01, 'ro'))

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
        
        # dummy_panel=np.zeros(len(panel))
        # for i in range(len(panel)):
        # 	dummy_panel[i]=panel[i].gamma1
        
        # print dummy_panel[0],panel[0].gamma1
        panel = solve_panel_method(flow_elements, panel)
        rhs_2 = compute_flowelement_velocity(dummy_flow_elements, panel)
        
        
        # for i in range(len(dummy_panel)):
        	# panel[i].gamma1=dummy_panel[i]
        
        # print dummy_panel[0],panel[0].gamma1
        
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

def compute_near_velocity(flow_elements, t):
    # location = np.zeros(nx*ny,dtype=complex)
    nx,ny = 40,40
    y,x = np.mgrid[-2:2:40j, -2:2:40j]
    # velocity = np.zeros((nx,ny),dtype=complex)
    U = np.zeros((nx,ny));
    V = np.zeros((nx,ny));
    for i in range(ny):
        for j in range(nx):
            location = complex(x[i,j],y[i,j])
            # print location
            temp = compute_total_velocity(body, flow_elements, panel, location)
            U[i,j] = temp.real
            V[i,j] = temp.imag


    print np.shape(U)
    # print np.asarray(velocity.real)
   	# U = velocity.real
   	# V = velocity.imag
    plt.plot(body.nodes.real, body.nodes.imag, linestyle='-')
    plt.streamplot(x, y, U, V, color=U, linewidth=2, cmap=plt.cm.autumn)
    plt.colorbar()
    # plt.quiver(location.real, location.imag, velocity.real, velocity.imag)
    
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlim(-1.5, 3)
    plt.ylim(-1.5, 1.5)
    filename = "./Velocity_field/%04d.png" %(t)
    plt.savefig(filename)
    plt.clf()
    # plt.show()

###########################################################################################
def solve_linear_system(flow_elements):
    b = np.zeros(body.n_panels+1)
    for i in range(body.n_panels):
        v = complex(0,0)
        v = compute_vortex_field(flow_elements, body.cp[i]) + free_stream.field()
        b[i] = -1*dot_product(v, body.cp[i])

    b[-1] = circulation
    
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

def init_panels(flow_elements):
    panel = []
    gamma = solve_linear_system(flow_elements)
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
	# total_velocity += image_vortex.blob_field(location)
    return total_velocity

###########################################################################################

def compute_flowelement_velocity(flow_elements, panel):
    flow_element_velocity = np.zeros(len(flow_elements),dtype=complex)
    for i in range(len(flow_elements)):
        location = flow_elements[i].current_location
        flow_element_velocity[i] = compute_total_velocity(body, flow_elements, panel, location) 
    return flow_element_velocity

###########################################################################################

def compute_slip(body, flow_elements, panel):
    slip_velocity = np.zeros(body.n_panels)
    for i in range(body.n_panels):
        total_velocity = complex(0,0)
        location = (1)*cmath.exp(body.cp_theta[i]*1j)
        # location = body.cp[i]
        total_velocity = compute_total_velocity(body, flow_elements, panel, location)
        slip_velocity[i] = dot_product(total_velocity,body.tangent[i])
    return slip_velocity

###########################################################################################

def distribute_vortex_blobs(body, flow_elements, slip_velocity):
    blob_radius = body.panel_length/np.pi
    blob_radial_pos = np.linalg.norm(1) + blob_radius
    image_blob_strength = 0
    for i in range(body.n_panels):
        new_blobs = []
        blob_strength_max = 0.1
        color = 'ro'
        blob_strength = slip_velocity[i]*body.panel_length
        image_blob_strength += blob_strength
        if blob_strength<0:
            blob_strength_max *= -1
            color = 'bo'
        blob_theta = body.cp_theta[i]
        blob_position = blob_radial_pos*cmath.exp(blob_theta*1j)
        blob_quant = int(abs(blob_strength/blob_strength_max))

        for j in range(blob_quant):
            new_blobs.append(Vortex(blob_strength_max, blob_position, blob_radius, color))
        remaining_strength = blob_strength - blob_quant*blob_strength_max
        new_blobs.append(Vortex(remaining_strength, blob_position, blob_radius, color))
        
        flow_elements += new_blobs
    # print image_blob_strength
###########################################################################################

def diffuse(flow_elements):
    current_positions = get_position_array(flow_elements)
    current_positions.real += np.random.randn(len(current_positions))*sigma + mean
    current_positions.imag += np.random.randn(len(current_positions))*sigma + mean
    current_positions = mirror(current_positions)
    store_positions(flow_elements, current_positions)

###########################################################################################

def mirror(current_positions):
    check_location = np.abs(current_positions)<np.abs(body.radius)
    particles_outside = list(compress(xrange(len(check_location)), check_location))
    
    if len(particles_outside)!=0:
        for item in particles_outside:
            old_radius = np.abs(current_positions[item])
            old_theta = cmath.phase(current_positions[item])
            
            new_radius = old_radius + 2*(body.radius-old_radius)
            new_theta = old_theta
            
            current_positions[item] = new_radius*cmath.exp(new_theta*1j)
    return current_positions

###########################################################################################

def plot_particles(flow_elements, t):
    plt.plot(body.nodes.real, body.nodes.imag,'k-')
    for i in range(len(flow_elements)):
        plt.plot(flow_elements[i].current_location.real, flow_elements[i].current_location.imag,flow_elements[i].color,markersize=2)

    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlim(-2, 6)
    plt.ylim(-2, 2)
    filename = "./Case_2/%04d.png" % (t)
    plt.savefig(filename)
    plt.clf()

###########################################################################################

def solve_panel_method(flow_elements, panel):
    gamma = solve_linear_system(flow_elements)
    panel = reset_panels(panel, gamma)
    return panel

###########################################################################################

def convect(flow_elements, panel):
    flow_element_velocity = compute_flowelement_velocity(flow_elements, panel)
    flow_element_position = get_position_array(flow_elements)
    # print "before:",panel[0].gamma1
    new_flow_element_position = sol.rk2(flow_elements, panel, flow_element_velocity, flow_element_position)
    # print "Afer:",panel[0].gamma1
    store_positions(flow_elements, new_flow_element_position)

###########################################################################################

def compute_vortex_momentum(flow_elements, vortex_momentum):
	particle_momentum = 0
	for vortex in flow_elements:
		particle_momentum += vortex.strength*complex(vortex.current_location.imag,-vortex.current_location.real)
	particle_momentum *= rho
	vortex_momentum.append(particle_momentum)

###########################################################################################

def compute_drag_force(vortex_momentum):
	return 0;

###########################################################################################
def solve_flow(flow_elements, panel):
    for t in range(nt):
        if t%20==0:
            print(str(int(100.0*t/nt))+"%........")

        #Solve Panel Method for no-penetration
        panel = solve_panel_method(flow_elements, panel)

        # total_gamma = image_vortex.strength
        # for item in panel:
        # 	total_gamma += item.gamma1
        # for item in flow_elements:
        # 	total_gamma += item.strength
        # print total_gamma

        

        #compute and store slip velocity at the control points
        slip_velocity = compute_slip(body, flow_elements, panel)
        # print slip_velocity

        #convect the flow elements
        convect(flow_elements, panel)
        
        #plot the quiver plot of velocity
        # compute_near_velocity(flow_elements, t)


        #distribute vortex blobs at the control points
        distribute_vortex_blobs(body, flow_elements, slip_velocity)
        print t, len(flow_elements)
        
        #Diffuse all the flow elements
        diffuse(flow_elements)

       	#Calculate vortex momentum
       	# compute_vortex_momentum(flow_elements, vortex_momentum)
       	# print vortex_momentum
        
        #plot the particles at every time step
        plot_particles(flow_elements, t)

        
        
    print "[DONE]"
    
###########################################################################################

#initiate Flow elements
flow_elements = []
vortex_momentum = []
free_stream = Uniform(v_infinity,0)
# image_vortex = Vortex(0,complex(0,0),0.04,'bo')

###########################################################################################

#initiate solver
sol = Solver(dt)

###########################################################################################

#initiate panels
panel = init_panels(flow_elements)


###########################################################################################
# test_single_vortex_motion(flow_elements)
#Solve for the time evolving system
solve_flow(flow_elements, panel)

###########################################################################################

def plot_tracers(location_tracer_time, option, file_name):
    if option=='movie':
        theta = np.zeros((4, 1))
        radius_vector = np.zeros((4, 2), dtype=complex)
        for i in range(len(theta)):
            theta[i] = i*np.pi/2
            
        for t in range(nt):
            if t%5==0:
                plt.plot(body.nodes.real, body.nodes.imag, linestyle='-')
                for l in range(len(theta)):
                    radius_vector[l, 0] = cmath.exp(1j*theta[l][0])
                    plt.plot(radius_vector[l, :].real, radius_vector[l, :].imag, linestyle='-')
                
                for j in range(len(location_tracer_time[0, :])):
                    if t!=0:
                        plt.plot(location_tracer_time[0:t-1, j].real, location_tracer_time[0:t-1, j].imag, linestyle='-')
                    else:
                        plt.plot(location_tracer_time[0, j].real, location_tracer_time[0, j].imag, linestyle='-')
                    plt.plot(location_tracer_time[t, j].real, location_tracer_time[t, j].imag, 'o')
                plt.gca().set_aspect('equal', adjustable='box')
                plt.xlim(-6, 6)
                plt.ylim(-6, 6)
                filename = "./Case_1/%04d.png" % (t/5)
                plt.savefig(filename)
                plt.clf()
            theta = theta + dt*omega

    elif option=='pathlines':
        plt.plot(body.nodes.real, body.nodes.imag, linestyle='-')
        for j in range(len(location_tracer_time[0, :])):
            plt.plot(location_tracer_time[:, j].real, location_tracer_time[:, j].imag, linestyle='-')
            plt.plot(location_tracer_time[-1, j].real, location_tracer_time[-1, j].imag, 'o')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlim(-6, 6)
        plt.ylim(-6, 6)
        plt.xlabel('X coordinate')
        plt.ylabel('Y coordinate')
        filename = "./"+file_name+".png"
        plt.savefig(filename)
        plt.clf()
    else:
        print "Wrong plotting option"
###########################################################################################

def plot_general(location_flowelement_time, filename):
    plt.plot(body.nodes.real, body.nodes.imag, linestyle='-')
    plt.plot(location_flowelement_time[:, 0].real, location_flowelement_time[:, 0].imag, linestyle='-')

    plt.plot(location_flowelement_time[-1, 0].real, location_flowelement_time[-1, 0].imag, 'o')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlim(-6, 6)
    plt.ylim(-6, 6)
    plt.savefig('./'+filename+'.png')
    plt.show()










