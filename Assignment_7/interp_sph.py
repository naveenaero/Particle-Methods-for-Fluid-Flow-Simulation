import numpy as np
import matplotlib.pyplot as plt
from math import sin,cos,pi,exp,sqrt

test_function = lambda x: np.sin(2*pi*x)
test_function_derivative = lambda x: 2*pi*np.cos(2*pi*x)
gaussian_kernel = lambda q,h: np.exp(-(q/h)**2)/(sqrt(pi)*h)
gaussian_kernel_derivative = lambda q,h: -2*np.sign(q)*(abs(q)/h**2)*gaussian_kernel(q,h)



def cubic_spline_kernel(rij, h):
    ''' Function to implement the cubic spline
        smoothing kernel
    '''
    r = abs(rij)
    q = r/h
    q2 = q**2
    q3 = q**3

    if q < 1.0 :
        func = (2./3 - q2 + 0.5*q3)/h
        return func
    elif q >= 1.0 and q < 2.0:
        func = (2-q)**3/(6*h)
        return func
    else:
        return 0

def cubic_spline(rij, h):
    ''' Function to implement the cubic spline
        smoothing kernel
    '''
    r = np.abs(rij)
    q = r/h
    q2 = q*q
    q3 = q*q*q

    cond1 = q<1.0
    cond2 = (q>=1.0) & (q<2.0)
    cond3 = q>=2.0

    func1 = (2./3 - q2 + 0.5*q3)/h
    func2 = (2-q)**3/(6*h)
    func3 = 0

    return cond1*func1 + cond2*func2 + cond3*func3 

    
def cubic_spline_derivative(rij, h):
    ''' Function to implement the cubic spline
        smoothing kernel derivative
    '''
    r = np.abs(rij)
    q = r/h
    r2 = r**2
    h2 = h*2
    h3 = h**3

    if q < 1.0:
        func = (-2*r/h2 + 1.5*r2/h3)/h
        return np.sign(rij)*func
    elif q >= 1.0 and q < 2.0:
        func = -1*(2-q)**2/(2*h2)
        return np.sign(rij)*func
    else:
        return 0

def cubic_spline_kernel_derivative(rij, h):
    ''' Function to implement the cubic spline
        smoothing kernel derivative
    '''
    r = np.abs(rij)
    q = r/h
    r2 = r*r
    h2 = h*h
    h3 = h*h*h
    
    cond1 = q<1.0
    cond2 = (q>=1.0) & (q<2.0)
    cond3 = q>=2.0

    func1 = (-2*r/h2 + 1.5*r2/h3)/h
    func2 =  -1 *(2-q)**2/(2*h2)
    func3 = 0

    return np.sign(rij)*(cond1*func1 + cond2*func2 + cond3*func3)

def test_spline():
    ''' Function to test the correct implementation
        of the spline smoothing kernel
    '''
    
    x = np.linspace(-3,3,100)
    y = np.zeros(len(x))
    y_der = np.zeros(len(x))
    h = 1
    
    for i in range(len(x)):
        y[i] = cubic_spline_kernel(x[i],h)
        y_der[i] = cubic_spline_derivative(x[i],h)
    
    plt.plot(x, y, label='kernel')
    plt.plot(x, y_der, label='kernel_derivative')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('q')
    plt.title('Cubic-Spline Kernel')
    plt.show()

def test_gaussian():
    ''' Function to test the correct implementation
        of the gaussian smoothing kernel
    '''
    x = np.linspace(-3,3,100)
    y = np.zeros(len(x))
    y_der = np.zeros(len(x))
    h = 1
    for i in range(len(x)):
        y[i] = gaussian_kernel(x[i],h)
        y_der[i] = gaussian_kernel_derivative(x[i],h)

    plt.plot(x, y, label='kernel')
    plt.plot(x, y_der, label='kernel_derivative')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('q')
    plt.title('Gaussian Kernel')
    plt.show()

def test_kernel():
    test_gaussian()
    test_spline()

def plot_error(L2_error_n, n, L2_error_h, h):
    plt.title('Error variation with number of base points and "h" value')
    plt.subplot(2,1,1)
    plt.plot((np.asarray(n)), np.log10(np.asarray(L2_error_n)), '-o')
    plt.xlabel('Number of base points')
    plt.ylabel('log(L2 error norm)')
    
    plt.subplot(2,1,2)
    plt.plot((np.asarray(h)), np.log10(np.asarray(L2_error_h)), '-o')
    plt.xlabel('h')
    plt.ylabel('log(L2 error norm)')
    plt.show()
    # plt.savefig('Error_plot_no_derivative')
    

def plot_function(x, y, x1, y1, kernel):
    ''' plot the approximated and exact functions
    ''' 
    # plt.subplot(2,1,1)
    plt.plot(x, y,  label='SPH-approximated-'+kernel)
    # plt.title("SPH approximated - " + kernel)
    # plt.ylabel('<f(x)>')
    # plt.xlabel('x')
    
    # plt.subplot(2,1,2)
    plt.plot(x1, y1, 'ro', label='exact')
    # plt.title('Exact function')
    plt.ylabel('f(x);  <f(x)>')
    plt.xlabel('x')
    plt.xlim(np.min(x)-0.5,np.max(x)+0.5)
    plt.ylim(np.min(y)-0.5,np.max(y)+0.5)
    plt.title(kernel+' SPH approximation')
    plt.legend()
    


def interp_cubic_spline_kernel(x_interp, x_base, y_base, h, kernel, derivative):
    ''' compute the SPH aproximation (x,y) of the 
        given function points at the data points (x1,y1)
        using the cubic spline kernel
    '''
    dx = x_base[1]-x_base[0]
    n_interp = len(x_interp)
    n_base = len(x_base)
    fx = np.zeros(n_interp)

    if kernel == 'cubic-spline':
        if derivative:
            for i in range(n_interp):
                for j in range(n_base):
                    fx[i] += dx*y_base[j]*cubic_spline_derivative(x_interp[i]-x_base[j], h)
            return fx
        else:
            for i in range(n_interp):
                for j in range(n_base):
                    fx[i] += dx*y_base[j]*cubic_spline_kernel(x_interp[i]-x_base[j], h)
            return fx
    elif kernel=='gaussian':
        if derivative:
            for i in range(n_interp):
                for j in range(n_base):
                    fx[i] += dx*y_base[j]*gaussian_kernel_derivative(x_interp[i]-x_base[j], h)
            return fx
        else:
            for i in range(n_interp):
                for j in range(n_base):
                    fx[i] += dx*y_base[j]*gaussian_kernel(x_interp[i]-x_base[j], h)
            return fx

def compute_error_norm(true_sol, computed_sol):
    ''' Function to compute the L2 error norm
    '''
    return np.linalg.norm(computed_sol-true_sol)/np.linalg.norm(true_sol)


def calculate_errors(kernel, derivative):
    L2_err_n = []
    L2_err_h = []
    n = []
    h = []
    xmin = -1.0
    xmax = 1.0
    if kernel=='cubic-spline' and derivative:
        h_factor = 1.0
    else:
        h_factor = 1.2
    for i in range(8):
        # change the number of base points
        n.append(10*2**i)
        # interpolate at 100 points for different base-points
        x_interp = np.linspace(-1,1,100)
        # calculate error
        L2_err_n.append(interp_sph(x_interp, n[-1], h_factor*(xmax-xmin)/(n[-1]-1), kernel, derivative))

    
    h_factor = [0.6, 0.7, 0.72, 0.75, 0.77, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]
    interp_points = 100
    base_points = 100
    x_interp = np.linspace(xmin,xmax,interp_points)
    
    for i in range(len(h_factor)):
        h.append(h_factor[i]*(xmax-xmin)/(base_points-1))
        L2_err_h.append(interp_sph(x_interp, base_points, h[-1], kernel, derivative))

    plot_error(L2_err_n, n, L2_err_h, h)



def interp_sph(x_interp, n_base, h, kernel, derivative):
    ''' Function to approximate a function using
        the SPH smoothing kernel
    '''
    #define interpolation points - equi-spaced - domain (-1,1)
    x_base = np.linspace(np.min(x_interp), np.max(x_interp), n_base)
    y_base = test_function(x_base) + np.random.randn(len(x_base))*sqrt(0.000)
    
    #interpolate using the mentioned smoothing kernel - Gaussian/cubic spline
    y_interp = interp_cubic_spline_kernel(x_interp, x_base, y_base, h, kernel, derivative)
    
    #compute L2 error norm
    if derivative:
        L2_err = compute_error_norm(test_function_derivative(x_interp), y_interp)
    else:
        L2_err = compute_error_norm(test_function(x_interp), y_interp)
    
    print "L2-error = ", L2_err
    
    #plot the results
    if derivative:
        plot_function(x_interp, y_interp, x_base, test_function_derivative(x_base), kernel)
    else:
        plot_function(x_interp, y_interp, x_base, y_base, kernel)
    plt.show()
    
    return L2_err


def test_interp(xmin, xmax, n_base, n_interp, h_factor):
    dx = float(xmax-xmin)/(n_base-1)
    h = h_factor*dx
    x_interp = np.linspace(xmin, xmax, n_interp)
    interp_sph(x_interp, n_base, h, 'gaussian', False)
    


