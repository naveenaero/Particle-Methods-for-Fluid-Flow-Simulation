import numpy as np
import math

number_of_vortex_elements = 100

gamma0 = 1
b = 1
Gamma = lambda y:gamma0*4*y/(b**2*math.sqrt(1-4*y**2/b**2))

spanwise_locations = np.linspace(-0.499,0.499,number_of_vortex_elements)
dx = spanwise_locations[1]-spanwise_locations[0]

f = open('flow_elements.txt','r+')
f.truncate()

#Write info in text file
for i in range(len(spanwise_locations)):
    f.write('vortex 0 '+str(spanwise_locations[i])+' '+str(dx*Gamma(spanwise_locations[i]))+'\n');

f.close()
     
