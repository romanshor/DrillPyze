from DrillPyze import well, drillstring, statics, dynamics, util
import numpy as np
import matplotlib.pyplot as plt
import time


# Set up parameters
simulation_length = 1200             # seconds
bit_depth = 7000                    # ft
casing_depth = 1000                 # ft
bit_size = 5                        # in


# Topdrive RPM setpoint commands as an matrix in the format (time (s), RPM)
TD_Setpoints = [[0,0], [1,30], [20,45], [50, 60], [360, 75]]
setpoint_idx = 0                    # helper variable - current index in setpoint matrix
TD_Setpoint = 0                     # helper variable - current rpm setpoint

# bitconstant changes in the format (time (s), constant)
rock_constants = [[0,0], [20, 1000], [50, 2500], [200, 3000], [250, 1000], [350, 1250], [400, 3500], [800, 2500]]
rock_idx = 0
rock_constant = 0

# Import the wellbore and BHA
well1 = well.wellbore({'xlsx':'../data/survey.xlsx'}, casing_depth, bit_size)
bha1 = drillstring.BHA('BHA 1', '../data/bha.xlsx')

# Set up the simulation
sim = dynamics.Simulation('Test Simulation', well1, bha1)
sim.initialize_simulation(bit_depth, simulation_length, bit_constant=1750)

print(f"Initializing Simulation with a simulation length of {simulation_length:0.1f} seconds.")

# Do the steadystate calculations
(Ft, Ft_pu, Ft_so, Fn, M) = statics.normal_forces(well1, bha1)

# Set up the figure
(fig, axes) = util.setup_plotting_estimates(live=True)

# define the number of frames per second
plot_frequency = 10

print(f"Iterating at {plot_frequency:0.1f} fps...")


data = np.genfromtxt('../data/data_long.csv', delimiter=',')
sim.initialize_estimator(len(data), intensity=250000)


tic = time.perf_counter()

# Loop the simulation
for i in range (0,simulation_length*plot_frequency):
    
    # check for the next TD Setpoint change
    if setpoint_idx < len(TD_Setpoints) and i == TD_Setpoints[setpoint_idx][0]*plot_frequency:
        TD_Setpoint = TD_Setpoints[setpoint_idx][1]
        setpoint_idx +=1


    # step the simulation forward by 1/plot_frequency
    sim.step(1/plot_frequency, TD_RPM_SP=TD_Setpoint)
    #sim = dynamics.Simulation.speedhelper(sim,1/plot_frequency, TD_RPM_SP=TD_Setpoint)

    sim.update_estimate(data[i,3])

# update the plots
util.plot_estimates(fig, axes, sim, bha1, M, data)



toc = time.perf_counter()
print(f"Estimation completed in {toc - tic:0.4f} seconds") 







