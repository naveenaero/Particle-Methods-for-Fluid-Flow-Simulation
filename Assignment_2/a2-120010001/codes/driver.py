from compute_flow import *

option  = input('Want the Flow Element to move(1) or not move(0)?:')
method = ['rk2']

time_horizon = 1
dt0 = 0.001
dt = np.array([[dt0]])#[dt0/2],[dt0/4],[dt0/8]])
	
error = np.zeros((len(dt),2))
num = 0


'''
code snippet to get the solution for various schemes and for different time step size
'''
for item in method:
	solution_method = item 

	for u in range(len(dt)):
		[TracerElements,Elements] = get_flow_elements()
		NumberOfTracers = len(TracerElements)
		NumberOfElements = len(Elements)

		# AnalayticalSolution = compute_analytical_solution(time_horizon,dt[u],Elements)

		nt = int(time_horizon/dt[u])

		sol = Solver(dt[u])
		updater = FlowUpdate(option)
		t = np.linspace(0,time_horizon,nt+1)


		ElementPosition = np.zeros((nt+1,NumberOfElements))*complex(0,0)
		ElementPosition = store_position_data(Elements,ElementPosition,0)


		TracerPosition = np.zeros((nt+1,NumberOfTracers))*complex(0,0)
		TracerPosition = store_tracer_data(TracerElements,TracerPosition,0)

		'''
		Loop for Time Integration
		'''
		for k in range(nt):
			##############################################################################
			# Solve for Tracer motion
			if NumberOfTracers>0:
				TracerVelocity = updater.get_tracer_velocity(TracerElements,Elements)
				curr_pos = convert_to_array(Elements,TracerElements,'Tracer')
				# pass instances of class Solver (sol) and FlowUpdate (updater) for caluclating new location
				# flag=1 computes for FlowElements, flag=0 computes for Tracers
				NewLocation = updater.update_position(sol,updater,TracerElements,Elements,curr_pos,TracerVelocity,solution_method,0)
				TracerElements = convert_from_array(Elements,TracerElements,NewLocation,'Tracer')
				TracerPosition = store_tracer_data(TracerElements,TracerPosition,k+1)
			
			##############################################################################
			# Solve for Flow ELement motion
			ElemVelocity = updater.get_flow_element_velocity(Elements)
			curr_pos = convert_to_array(Elements,TracerElements,'Elements')
			NewLocation = updater.update_position(sol,updater,TracerElements,Elements,curr_pos,ElemVelocity,solution_method,1)
			Elements = convert_from_array(Elements,TracerElements,NewLocation,'Elements')
			ElementPosition = store_position_data(Elements,ElementPosition,k+1)
			
			##############################################################################

			
		# calculate the relative error
		# error[u,num] = np.linalg.norm(abs(AnalayticalSolution[-1]-ElementPosition[-1,1]))/np.linalg.norm(abs(AnalayticalSolution[-1]))

	num = num+1



'''
code snippet for plotting the error convergence plots
'''
################################################################################
# fig, ax = plt.subplots()
# ax.loglog(dt,error[:,0],label='euler')
# ax.loglog(dt,error[:,1],label='rk2')
# legend = ax.legend(loc='upper left')

# # The frame is matplotlib.patches.Rectangle instance surrounding the legend.
# frame = legend.get_frame()
# frame.set_facecolor('0.90')

# # Set the fontsize
# for label in legend.get_texts():
#     label.set_fontsize('large')

# for label in legend.get_lines():
#     label.set_linewidth(1.5)  # the legend line width
# plt.xlabel('Log(dt)')
# plt.ylabel('Log(error)')
# slope_euler,intercept_euler = np.polyfit(dt,error[:,0],1)
# slope_rk2,intercept_rk2 = np.polyfit(dt,error[:,1],1)
# print ("euler_slope="),slope_euler
# print ("rk2_slope="),slope_rk2
# plt.savefig('convergence.png')
###############################################################################



'''
code snippet for plotting Flow element position history
'''
##############################################################################
# for j in range(NumberOfElements):
# 	plt.plot(ElementPosition[:,j].real,ElementPosition[:,j].imag)

plt.plot(ElementPosition[-1,:].real,ElementPosition[-1,:].imag,linestyle='-')
plt.xlim(-1, 1)
plt.ylim(-1, 1)
plt.xlabel('X coordinate')
plt.ylabel('Y coordinate')
plt.gca().set_aspect('equal', adjustable='box')
plt.title('Vortex sheet roll-up \t Time Horizon='+str(time_horizon)+' sec')
plt.savefig("./with_krasney_blob.png")
plt.show()
plt.clf()
##############################################################################



'''
code snippet for plotting tracer element position history
'''
##############################################################################
# for j in range(NumberOfTracers):
# 	plt.plot(TracerPosition[:,j].real,TracerPosition[:,j].imag)
# 	plt.plot(TracerPosition[-1,j].real,TracerPosition[-1,j].imag,'ro')
# for j in range(NumberOfElements):
# 	if Elements[j].kind!='uniform':
# 		plt.plot(ElementPosition[0,j].real,ElementPosition[0,j].imag,'ro')
# plt.xlim(-4, 4)
# plt.ylim(-2, 2)
# plt.xlabel('X coordinate')
# plt.ylabel('Y coordinate')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.title('Tracer + Flow Elements')
# #plt.savefig('source_sink_tracer_new.png')
# plt.show()
# plt.clf()
##############################################################################







