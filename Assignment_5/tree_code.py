import numpy as np
import matplotlib.pyplot as plt
from itertools import compress
import matplotlib.patches as patches
import time
import math



def get_random_vortices(n_vortices):
    pos = np.zeros(n_vortices, dtype=complex)
    pos.real = np.random.random(n_vortices)
    pos.imag = np.random.random(n_vortices)
    strength = np.random.random(n_vortices)
    # pos.real = np.random.normal(0,0.4,n_vortices)
    # pos.imag = np.random.normal(0,0.4,n_vortices)
    # strength = np.random.normal(1,0.5,n_vortices)
    blob_radius = 1e-4
    flow_elements = []
    for i in range(n_vortices):
        flow_elements.append(cf.Vortex(strength[i], pos[i], blob_radius, 'bo'))
    return flow_elements



def get_position_array(flow_elements):
    curr_pos = []
    for i in range(len(flow_elements)):
        curr_pos.append(flow_elements[i].current_location)
    return np.asarray(curr_pos)
 

def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)

def n_squared(flow_elements):
    n_vortices = len(flow_elements)
    velocity = np.zeros(n_vortices, dtype=complex)
    for i in range(n_vortices):
        velocity[i] = cf.compute_vortex_field(flow_elements, flow_elements[i].current_location)
    return velocity

class Cell:
    def __init__(self, x, y, n, indices):
        # geometrix parameters of cell
        self.xmin = x[0]
        self.xmax = x[1]
        self.ymin = y[0]
        self.ymax = y[1]
        self.center = complex(np.average(x), np.average(y))
        self.radius = abs(self.center - complex(np.average(x), self.ymin))

        # number of vortices in cell
        self.n = n
        # indices of vortices inside the cell
        self.indices = indices
        # level of cell
        self.level = 0
        # index of cell
        self.index = 0
        # parent of cell
        self.parent = None
        # childs of cell
        self.children = None

        # order of multipole expansion
        self.order = 15
        # multipole coefficients
        self.a = np.zeros(self.order+1, dtype=complex)
        
    def create_children(self, vortices):
        '''create children of cell given all vortex positions'''
        self.children = []
        x = np.linspace(self.xmin, self.xmax, 3)
        y = np.linspace(self.ymin, self.ymax, 3)
        for i in range(2):
            for j in range(2):
                x_new = [x[i], x[i+1]]
                y_new = [y[j], y[j+1]]
                [n_inside, indices] = get_number_of_vortices(vortices, x_new, y_new)
                new_cell = Cell(x_new, y_new, n_inside, indices)
                new_cell.level = self.level + 1
                new_cell.parent = self
                self.children.append(new_cell)

        return self.children

    def compute_multipole_coefficients(self, flow_elements):
        '''compute multipole expansion of vortices at the cell center'''
        z_o = self.center
        for j in range(1,self.order+1):
            for k in self.indices:
                z_k = flow_elements[k].current_location
                rel_position = z_k - z_o
                gamma_k = flow_elements[k].strength
                self.a[j] += gamma_k * rel_position**(j-1)

    def transfer_multipole_coefficients(self):
        '''transfer multipole coefficients to parent'''
        if self.parent!=None:
            z_c = self.parent.center
            z_o = self.center
            z_co = z_o - z_c
            b = np.zeros(self.order+1, dtype=complex)
            for j in range(1, self.order+1):
                for k in range(1, j+1):
                    b[j] += self.a[k]*nCr(j-1, k-1)*(z_co)**(j-k)
            self.parent.a += b

    def multipole_expansion(self, z_p):
        '''compute multipole expansion at particle location z_p'''
        velocity = complex(0,0)
        z_relative = (z_p - self.center)
        for j in range(len(self.a)):
            velocity += self.a[j]/(z_relative)**j
        velocity *= complex(0, -1./(2*np.pi))
        return np.conj(velocity)

    def direct_velocity(self, flow_elements, z_p):
        '''compute direct velocity using order N-squared algorithm'''
        velocity = complex(0,0)
        for i in self.indices:
            velocity += flow_elements[i].blob_field(z_p)
        return velocity


def get_number_of_vortices(vortices, x, y):
    ''' get number of vortices and their global indices in a given boundary '''
    check_location = (vortices.real<x[1]) & (vortices.real>=x[0]) & (vortices.imag<y[1]) & (vortices.imag>=y[0])
    particles_inside = list(compress(xrange(len(check_location)), check_location))
    return len(particles_inside), particles_inside

def create_root(vortices):
    ''' create the root of the tree '''
    # get the min and max x & y locations of vortices
    d = 0.01
    x = [np.min(vortices.real)-d, np.max(vortices.real)+d]
    y = [np.min(vortices.imag)-d, np.max(vortices.imag)+d]
    delta_x = x[1]-x[0]
    delta_y = y[1]-y[0]
    width = max(delta_x, delta_y)
    x = [x[0], x[0]+width]
    y = [y[0], y[0]+width]
    # get the number of vortices inside along with their indices
    [n_inside, indices] = get_number_of_vortices(vortices, x, y)
    # create root
    root = Cell(x, y, n_inside, indices)
    # assign root its level and index
    root.level = 0
    root.index = 0
    return root

def create_tree(flow_elements, max_vortices):
    ''' generate the tree '''
    vortices = get_position_array(flow_elements)
    root = create_root(vortices)
    stack = [root]
    cells = [root]
    index = 1
    while len(stack)>0:
        # get the last cell in the stack
        cell = stack.pop()
        # compare number of vortices with max_vortices
        if cell.n > max_vortices:
            # if greater number of vortices found, divide cell into children
            children = cell.create_children(vortices)
            # assign level and index to childre
            for child in children:
                child.index = index
                index += 1
            # extend the stack and cell list with these newly created children
            stack.extend(children)
            cells.extend(children)
    # return cells
    return cells, vortices

def compute_multipoles(cells, n_levels, flow_elements):
    ''' compute multipoles of all cells '''
    # for all childless cells, compute multipole coefficients
    # cell bunch - cell list at a particular level
    for level in cells:
        # for cells at a certain level
        for cell in level:
            # check if it is childless
            if cell.children==None:
                # compute multipole coefficients of cell
                cell.compute_multipole_coefficients(flow_elements)
    
    # iterate on levels
    for level in range(n_levels-2, 1, -1):
        # iterate over cells at this level
        for cell in cells[level]:
            # iterate over children of a cell at this level if children exist
            if cell.children!=None:
                for child in cell.children:
                    # transfer multipole coefficients of child to its parent
                    child.transfer_multipole_coefficients()
           
def convert(cells, n_levels):
    ''' classify cells on basis of levels'''
    new_cells = []
    for i in range(n_levels):
        new_cells.append([])

    for cell in cells:
        new_cells[cell.level].append(cell)
    
    return new_cells

def get_number_of_levels(cells):
    n_levels = 0
    for cell in cells:
        if cell.level>n_levels:
            n_levels = cell.level
    return n_levels + 1

def compute_cell_velocity(cell, z_p, flow_elements):
    velocity = complex(0,0)
    if abs(cell.center - z_p) >= 2.0*cell.radius:
        velocity = cell.multipole_expansion(z_p)
    elif cell.children!=None:
        for child in cell.children:
            velocity += compute_cell_velocity(child, z_p, flow_elements)
    else:
        velocity = cell.direct_velocity(flow_elements, z_p)
    return velocity

def get_childless_cells(cells):
    childless = []
    for level in cells:
        for cell in level:
            if cell.children==None:
                childless.append(cell)
    return childless


def compute_velocity(cells, flow_elements):
    childless_cells = get_childless_cells(cells)
    velocity = np.zeros(len(flow_elements), dtype=complex)
    for cell in childless_cells:
        for i in cell.indices:
            p = flow_elements[i].current_location
            for b in cells[2]:
                velocity[i] += compute_cell_velocity(b, p, flow_elements)
            
    return velocity

def plot_cells(cells, vortices):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    for p in [
                patches.Rectangle(
                            (cells[i].xmin, cells[i].ymin),
                            cells[i].xmax-cells[i].xmin,
                            cells[i].ymax-cells[i].ymin,
                            fill=False
                            ) for i in range(len(cells))
                ]:
                    ax.add_patch(p)
    
    plt.plot(vortices.real, vortices.imag, 'bo', markersize=1)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Body Cell Treecode')
    plt.savefig('tree_code_cells.png')
    plt.clf()
    # plt.show()
    


def n_logn(flow_elements):
    # set maximum number of vortices allowed in a cell
    max_vortices = 20
    
    # create the tree
    [cells, vortices] = create_tree(flow_elements, max_vortices)
    
    # plot vortices and cells
    # plot_cells(cells, vortices)

    
    # compute the number of levels in the tree
    n_levels = get_number_of_levels(cells)
    
    # classify the cells on the basis of their level
    cells = convert(cells, n_levels)
    
    # compute the multipoles of the cells at their cell centres
    compute_multipoles(cells, n_levels, flow_elements)
    
    # compute the velocity of each particle using Tree code
    velocity = compute_velocity(cells, flow_elements)
    
    return velocity

def test(n_vortices):
    '''get both velocities'''
    flow_elements = get_random_vortices(n_vortices)
    
    t_tree=0; t_nsquared=0; t0=0; t1=0

    t0 = time.time()
    velocity_tree = n_logn(flow_elements)
    t1 = time.time()
    t_tree = float(t1-t0)
    
    t0 = time.time()
    velocity_nsquared = n_squared(flow_elements)
    t1 = time.time()
    t_nsquared = float(t1-t0)

    # print L2-error norm
    print "-------------------------------------"
    print "Number of Particles: ",n_vortices
    print "L2-error norm: ",np.linalg.norm(velocity_tree-velocity_nsquared)/np.linalg.norm(velocity_nsquared)
    print "Time for tree code: ", t_tree
    print "Time for N-squared: ", t_nsquared
    return float(t_tree), float(t_nsquared)

def time_plot():
    n_vortices = np.array([50,100,200,400,800])
    time_tree = np.zeros(len(n_vortices))
    time_squared = np.zeros(len(n_vortices))
    for n in range(len(n_vortices)):
        [a,b] = test(n_vortices[n])
        time_tree[n] = a
        time_squared[n]= b
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.plot(n_vortices, time_tree, '-o', label='Tree-code')
    plt.plot(n_vortices, time_squared, '-o', label='O(N^2)')
    plt.title('Comparison of Tree Code with O(N^2) algorithm')
    plt.xlabel('Number of Particles (N)')
    plt.ylabel('Time(s)')
    plt.legend()
    plt.savefig('comparison.png')
    # plt.show()

# time_plot()

