from DrillPyze import well, drillstring, statics, dynamics, util
import numpy as np
import matplotlib.pyplot as plt
import time

# Set up parameters
simulation_length = 240             # seconds
bit_depth = 7000                    # ft
casing_depth = 1000                 # ft
bit_size = 5                        # in

tic = time.perf_counter()

# Flags to save data and/or plot results while running
SAVE_DATA = False
PLOT_LIVE = False

# Topdrive RPM setpoint commands as an matrix in the format (time (s), RPM)
TD_Setpoints = [[0,0], [1,30], [20,45], [50, 60], [360, 75]]
setpoint_idx = 0                    # helper variable - current index in setpoint matrix
TD_Setpoint = 0                     # helper variable - current rpm setpoint

# bitconstant changes in the format (time (s), constant)
rock_constants = [[0,0], [10, 1000], [50, 2500], [200, 3000], [250, 1000], [350, 1250], [400, 3500], [800, 2500]]
rock_idx = 0                        # helper variable - ccurent index in rock_constant matrix
rock_constant = 0                   # helper variable - current bit-rock constant

# Import the wellbore and BHA
well1 = well.wellbore({'xlsx':'../data/survey.xlsx'}, casing_depth, bit_size)
bha1 = drillstring.BHA('BHA 1', '../data/bha.xlsx')

# Set up the simulation
sim = dynamics.Simulation('Test Simulation', well1, bha1)
sim.initialize_simulation(bit_depth, simulation_length)

print(f"Initializing Simulation with a simulation length of {simulation_length:0.1f} seconds.")

# Do the steadystate calculations
(Ft, Ft_pu, Ft_so, Fn, M) = statics.normal_forces(well1, bha1)

# Set up the figure
(fig, axes) = util.setup_plotting(live=PLOT_LIVE)

# define the number of frames per second
plot_frequency = 10
if PLOT_LIVE:
    print(f"Iterating at {plot_frequency:0.1f} fps...")
else:
    print(f"Iterating at with a time step of {sim.dt:0.5f} seconds...")

# set up the data array for saving data
if SAVE_DATA:
    data = np.zeros((simulation_length*plot_frequency,7))

# Loop the simulation
for i in range (0,simulation_length*plot_frequency):
    
    # check for the next TD Setpoint change
    if setpoint_idx < len(TD_Setpoints) and i == TD_Setpoints[setpoint_idx][0]*plot_frequency:
        TD_Setpoint = TD_Setpoints[setpoint_idx][1]
        setpoint_idx +=1

    # vary the bit-rock parameter
    if rock_idx < len(rock_constants) and i == rock_constants[rock_idx][0]*plot_frequency:
        rock_constant = rock_constants[rock_idx][1]
        rock_idx +=1


    # step the simulation forward by 1/plot_frequency
    sim.step(1/plot_frequency, TD_RPM_SP=TD_Setpoint, bit_constant=rock_constant)

    # update the plots
    if PLOT_LIVE:
        util.plot(fig, axes, sim, bha1, M)

    # save a line of data
    if SAVE_DATA:
        data[i,:] = [sim.time[sim.idx-1], TD_Setpoint, sim.omega_td[sim.idx-1], sim.tau_td[sim.idx-1], sim.omega_bit[sim.idx-1], sim.tau_bit[sim.idx-1], rock_constant]

# plot the final results
util.plot(fig, axes, sim, bha1, M)


# write the data to data.csv
if SAVE_DATA:
    np.savetxt('data_long.csv', data, delimiter=',')

# print the simulation time
toc = time.perf_counter()
print(f"Simulation completed in {toc - tic:0.4f} seconds") 
