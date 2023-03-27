"""
A module containing functions for plotting the results of the simulation and other
useful functions.
"""



import matplotlib.pyplot as plt
import numpy as np


def setup_plotting(live=False):
    """
    Setup the plotting

    Parameters
    ==========
    live : bool
        if True, then live animation is available
    """
    plt.ion()
    fig = plt.figure()
    ax1 = fig.add_subplot(3,4,(3,4))
    ax2 = fig.add_subplot(3,4,(7,8))
    ax3 = fig.add_subplot(3,4,(1,5))
    ax4 = fig.add_subplot(3,4,(2,6))


    return fig, (ax1, ax2, ax3, ax4)

def setup_plotting_estimates(live=False):
    """
    Setup the plotting

    Parameters
    ==========
    live : bool
        if True, then live animation is available
    """
    plt.ion()
    fig = plt.figure()
    ax1 = fig.add_subplot(3,4,(3,4))
    ax2 = fig.add_subplot(3,4,(7,8))
    ax3 = fig.add_subplot(3,4,(1,5))
    ax4 = fig.add_subplot(3,4,(2,6))
    ax5 = fig.add_subplot(3,4,(9,10))
    ax6 = fig.add_subplot(3,4,(11,12))
    


    return fig, (ax1, ax2, ax3, ax4, ax5, ax6)

def plot(fig, axes, sim, bha, M):
    """
    Update the plot

    Parameters
    ==========
    fig : matplotlib.pyplot.figure
        Figure to plot in
    axes : tuple of axes objects
    sim : DrillPyze.dynamics.simulation
        Simulation object
    bha : DrillPyze.drillstring.BHA
        BHA object
    M : numpy.array
        Array of static friction torques
    """

    (ax1, ax2, ax3, ax4) = axes
    ax1.clear()
    ax1.plot(sim.time, sim.omega_bit/np.pi*30)
    ax1.plot(sim.time, sim.omega_td/np.pi*30, 'b')
    ax1.plot(sim.time, sim.omega_td_sp/np.pi*30, 'k')
    ax1.set_ylim(-30,120)

    ax2.clear()
    ax2.plot(sim.time, sim.tau_td, 'b')
    ax2.plot(sim.time, sim.tau_bit, 'r')

    ax3.clear()
    ax3.plot(sim.tau_v[1:-1,2],bha.x,'b')
    ax3.plot(M+sim.tau_v[-1,1],bha.x, 'k')

    ax4.clear()
    #ax4.plot(sim.omega[:,2]/(np.pi)*60,bha1.x, 'b')
    ax4.plot(sim.omega_v[1:-1,1]/np.pi*30,bha.x,'b')
    ax4.set_xlim(-30,120)

    fig.canvas.draw()
    fig.show()
    plt.pause(0.01)

    fig.show()



def plot_estimates(fig, axes, sim, bha, M, data):

    """
    Update the plot

    Parameters
    ==========
    fig : matplotlib.pyplot.figure
        Figure to plot in
    axes : tuple of axes objects
    sim : DrillPyze.dynamics.simulation
        Simulation object
    bha : DrillPyze.drillstring.BHA
        BHA object
    M : numpy.array
        Array of static friction torques
    """

    (ax1, ax2, ax3, ax4, ax5, ax6) = axes
    ax1.clear()
    ax1.plot(data[:,0], data[:,1], 'g--')
    #ax1.plot(data[:,0], data[:,2]/np.pi*30, 'c--')
    ax1.plot(sim.time, sim.omega_bit/np.pi*30)
    ax1.plot(sim.time, sim.omega_td/np.pi*30, 'b')
    ax1.set_ylim(-30,120)
    ax1.set_title('RPM vs time')

    ax2.clear()
    ax2.plot(sim.time, sim.tau_td, 'b')
    ax2.plot(sim.time, sim.tau_bit, 'r')
    ax2.plot(data[:,0], data[:,3], 'c--')
    ax2.set_title('Torque vs time)')

    ax3.clear()
    ax3.plot(sim.tau_v[1:-1,2],bha.x,'b')
    ax3.plot(M+sim.tau_v[-1,1],bha.x, 'k')
    ax3.set_title('Torque vs depth')

    ax4.clear()
    #ax4.plot(sim.omega[:,2]/(np.pi)*60,bha1.x, 'b')
    ax4.plot(sim.omega_v[1:-1,1]/np.pi*30,bha.x,'b')
    ax4.set_xlim(-30,120)
    ax4.set_title('RPM vs depth')

    ax5.clear()
    ax5.plot(data[:,0], data[:,6], 'c--')
    ax5.plot(data[:,0], sim.bit_estimates[:,1], 'b')
    ax5.set_title('Bit Estimate vs time')

    ax6.clear()
    ax6.plot(data[:,0], data[:,4]/np.pi*30, 'c--')
    ax6.plot(sim.time, sim.omega_bit/np.pi*30, 'b')
    ax6.set_title('Bottomhole RPM vs time')


    fig.canvas.draw()
    fig.show()
    plt.pause(0.01)

    fig.show()


def plot_estimates_multiple(fig, axes, sims, bha, M, data):

    """
    Update the plot

    Parameters
    ==========
    fig : matplotlib.pyplot.figure
        Figure to plot in
    axes : tuple of axes objects
    sim : DrillPyze.dynamics.simulation
        Simulation object
    bha : DrillPyze.drillstring.BHA
        BHA object
    M : numpy.array
        Array of static friction torques
    """

    (ax1, ax2, ax3, ax4, ax5, ax6) = axes
    
    for ax in axes:
        ax.clear()

    for sim in sims:
        ax1.plot(data[:,0], data[:,1], 'g--')
        #ax1.plot(data[:,0], data[:,2]/np.pi*30, 'c--')
        ax1.plot(sim.time, sim.omega_bit/np.pi*30)
        ax1.plot(sim.time, sim.omega_td/np.pi*30, 'b')
        ax1.set_ylim(-30,120)
        ax1.set_title('RPM vs time')

        ax2.plot(sim.time, sim.tau_td, 'b')
        ax2.plot(sim.time, sim.tau_bit, 'r')
        ax2.plot(data[:,0], data[:,3], 'c--')
        ax2.set_title('Torque vs time)')

        ax3.plot(sim.tau_v[1:-1,2],bha.x,'b')
        ax3.plot(M+sim.tau_v[-1,1],bha.x, 'k')
        ax3.set_title('Torque vs depth')

        #ax4.plot(sim.omega[:,2]/(np.pi)*60,bha1.x, 'b')
        ax4.plot(sim.omega_v[1:-1,1]/np.pi*30,bha.x,'b')
        ax4.set_xlim(-30,120)
        ax4.set_title('RPM vs depth')

        ax5.plot(data[:,0], data[:,6], 'c--')
        ax5.plot(sim.bit_estimates[:,0], sim.bit_estimates[:,1], 'b.')
        ax5.set_title('Bit Estimate vs time')

        ax6.plot(data[:,0], data[:,4]/np.pi*30, 'c--')
        ax6.plot(sim.time, sim.omega_bit/np.pi*30, 'b')
        ax6.set_title('Bottomhole RPM vs time')


    fig.canvas.draw()
    fig.show()
    plt.pause(0.01)

    fig.show()


            