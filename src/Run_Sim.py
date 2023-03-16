from DrillPyze import well, drillstring, statics, dynamics
import numpy as np
import matplotlib.pyplot as plt



# Import the wellbore and BHA
well1 = well.wellbore({'xlsx':'../data/survey.xlsx'}, 1000, 5)
bha1 = drillstring.BHA('BHA 1', '../data/bha.xlsx')

# Set up the simulation
sim = dynamics.Simulation('Test Simulation', well1, bha1)
sim.initialize_simulation(7000, 60)

# Do the steadystate calculations
(Ft, Ft_pu, Ft_so, Fn, M) = statics.normal_forces(well1, bha1)

# Set up the figure
plt.ion()
fig = plt.figure()
ax1 = fig.add_subplot(2,4,(3,4))
ax2 = fig.add_subplot(2,4,(7,8))
ax3 = fig.add_subplot(2,4,(1,5))
ax4 = fig.add_subplot(2,4,(2,6))

plot_frequency = 10
vectorized = True

if vectorized:

    for i in range (0,60*plot_frequency):
        if i > 30*plot_frequency:
            sim.step_vectorized(1/plot_frequency, 60)
        elif i > 1:
            sim.step_vectorized(1/plot_frequency, 30)
        else:
            sim.step_vectorized(1/plot_frequency)

        
        ax1.clear()
        ax1.plot(sim.time, sim.omega_bit)
        ax1.plot(sim.time, sim.omega_td, 'b')
        ax1.plot(sim.time, sim.omega_td_sp, 'k')

        ax2.clear()
        ax2.plot(sim.time, sim.tau_td, 'b')
        ax2.plot(sim.time, sim.tau_bit, 'r')

        ax3.clear()
        ax3.plot(sim.tau_v[1:-1,2],bha1.x,'b')
        ax3.plot(M+sim.tau_v[-1,1],bha1.x, 'k')

        ax4.clear()
        #ax4.plot(sim.omega[:,2]/(np.pi)*60,bha1.x, 'b')
        ax4.plot(sim.omega_v[1:-1,1],bha1.x,'b')

        fig.canvas.draw()
        fig.show()
        plt.pause(0.01)

    fig.show()

else:

    for i in range (0,60*plot_frequency):
        if i > 30*plot_frequency:
            sim.step(1/plot_frequency, 60)
        elif i > 1:
            sim.step(1/plot_frequency, 30)
        else:
            sim.step(1/plot_frequency)

        
        ax1.clear()
        ax1.plot(sim.time, sim.omega_bit)
        ax1.plot(sim.time, sim.omega_td, 'b')
        ax1.plot(sim.time, sim.omega_td_sp, 'k')

        ax2.clear()
        ax2.plot(sim.time, sim.tau_td, 'b')
        ax2.plot(sim.time, sim.tau_bit, 'r')

        ax3.clear()
        ax3.plot(sim.tau[:,2],bha1.x,'b')
        ax3.plot(M+sim.tau[-1,1],bha1.x, 'k')

        ax4.clear()
        #ax4.plot(sim.omega[:,2]/(np.pi)*60,bha1.x, 'b')
        ax4.plot(sim.omega[:,1],bha1.x,'b')

        fig.canvas.draw()
        fig.show()
        plt.pause(0.01)

    fig.show()

