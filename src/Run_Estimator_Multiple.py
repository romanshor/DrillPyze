from DrillPyze import well, drillstring, statics, dynamics, util
import numpy as np
import matplotlib.pyplot as plt
import time


# Set up parameters
simulation_length = 60             # seconds
bit_depth = 7000                    # ft
casing_depth = 1000                 # ft
bit_size = 5                        # in


# Import the wellbore and BHA
well1 = well.wellbore({'xlsx':'../data/survey.xlsx'}, casing_depth, bit_size)
bha1 = drillstring.BHA('BHA 1', '../data/bha.xlsx')

# Set up the simulations

sim1 = dynamics.Simulation('Test Simulation', well1, bha1)
sim1.initialize_simulation(bit_depth, simulation_length, bit_constant=1000)

sim2 = dynamics.Simulation('Test Simulation', well1, bha1)
sim2.initialize_simulation(bit_depth, simulation_length, bit_constant=2000)

sim3 = dynamics.Simulation('Test Simulation', well1, bha1)
sim3.initialize_simulation(bit_depth, simulation_length, bit_constant=3000)

sims = (sim1, sim2, sim3)
print(f"Initializing Simulation with a simulation length of {simulation_length:0.1f} seconds.")

# Do the steadystate calculations
(Ft, Ft_pu, Ft_so, Fn, M) = statics.normal_forces(well1, bha1)

# Set up the figure
(fig, axes) = util.setup_plotting_estimates(live=True)

# define the number of frames per second
plot_frequency = 10

print(f"Iterating at {plot_frequency:0.1f} fps...")


data = np.genfromtxt('../data/data_long.csv', delimiter=',')
datarate = 10
for sim in sims:
    sim.initialize_estimator(len(data), intensity=2500*plot_frequency)


tic = time.perf_counter()

# Loop the simulation
for i in range (0,simulation_length*plot_frequency):
    


    for sim in sims:
        # step the simulation forward by 1/plot_frequency
        sim.step(1/plot_frequency, TD_RPM_SP=data[int(i*datarate/plot_frequency),2])
        #sim = dynamics.Simulation.speedhelper(sim,1/plot_frequency, TD_RPM_SP=TD_Setpoint)

        sim.update_estimate(data[int(i*datarate/plot_frequency),3])

# update the plots
util.plot_estimates_multiple(fig, axes, sims, bha1, M, data)



toc = time.perf_counter()
print(f"Estimation completed in {toc - tic:0.4f} seconds") 







