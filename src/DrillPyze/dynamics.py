"""
A module that contains the simulation class and associated functions to enable
simulation of a coupled axial-torsional drillstring model with a choice of bit-rock
interaction models.

The estimation class is used to estimate the bit-rock interaction model parameters.
"""

import pandas as pd
import numpy as np
from DrillPyze import *
import matplotlib.pyplot as plt
from drawnow import drawnow
from pathos.multiprocessing import ThreadingPool as Pool  #@TODO: get ProcessingPool to work (10x speedup possible)
from scipy.special import expit
from numba import jit


class Simulation:
    """
    A class to maintain simulation state and describe a coupled axial-torsional 
    drillstring model with a choice of bit-rock interaction models
    """

    def __init__(self, name, wellbore, bha):
        """
        Create a new simulation with a bha and a wellbore
        
        Parameters
        ==========
        name : str
            Name of the simulation
        wellbore : DrillPyze.well.wellbore
            A wellbore object initialized with a well profile.
        hba : DrillPyze.drillstring.BHA
            A BHA object initialized with a drillstring configuration.

        """
        self.name = name
        self.wellbore = wellbore
        self.bha = bha
        self.dx = 90
        self.dt = 0.0001
        self.error = 0


    def initialize_simulation(self, bitdepth, simulation_length=100, update_period=30, bit_constant = 2000, params=None):
        """
        Initialize the simulation with the bit at a certain bitdepth
        
        Parameters
        ==========
        bitdepth : float
            Current depth of the bit
        simulation_length : int
            Length (in seconds) of the simulation.  Used to preallocate arrays.
            Default value is 100 seconds.        
        params : dict
            A dictionary of optional parameters to set up the simulation
            Default value is None
        """

        # Store the bit depth
        self.bitdepth = bitdepth
        

        # Discretize the Wellpath and Drillstring
        self.bha.discretize(self.dx, bitdepth)
        self.wellbore.get_interpolated(self.bha.x)


        # model parameters
        self.omega_c = 0.1                              # cut-off omega for static friction
        self.mu_s = 0.45                                # along-string static friction
        self.mu_k = 0.3                                 # along-string kinematic friction
        self.kt = 3.5                                   # inertial torsional damping
        self.ks = 0.3                                   # viscous damping coefficient

        # Do the static calculations on the drillstring
        (Ft, Ft_pu, Ft_so, Fn, M) = statics.normal_forces(self.wellbore, self.bha)
        self.Ft = Ft
        self.Fn = Fn
        self.M = M

        # Initialize storage for simulation variables

        # torsional model
        self.tau = np.zeros((len(self.bha.x),3))        # Torque along the string for t=t, t=t-dt and t=t-2dt
        self.omega = np.zeros((len(self.bha.x),3))      # RPM along the string for t=t, t=t-dt and t=t-2dt
        self.source = np.zeros((len(self.bha.x),1))     # source term along the drillstring

        self.tau_v = np.zeros((len(self.bha.x)+2,3))        # Torque along the string for t=t, t=t-dt and t=t-2dt.  Vectorized format
        self.omega_v = np.zeros((len(self.bha.x)+2,3))      # RPM along the string for t=t, t=t-dt and t=t-2dt.  Vectorized format
        self.source_v = np.zeros((len(self.bha.x)+2,1))     # Source term along the drillstring.  Vectorized format

        self.omega_td = np.zeros(int(simulation_length/self.dt))            # Top Drive RPM
        self.omega_td_sp = np.zeros(int(simulation_length/self.dt))         # Top Drive RPM Setpoint
        self.omega_bit = np.zeros(int(simulation_length/self.dt))           # Bit RPM
        self.tau_td = np.zeros(int(simulation_length/self.dt))              # Top Drive Torque
        self.tau_bit = np.zeros(int(simulation_length/self.dt))             # Bit Torque
        

        # axial model
        self.F = np.zeros((len(self.bha.x),3))          # Force along the string for t=t, t=t-dt and t=t-2dt
        self.v = np.zeros((len(self.bha.x),3))          # Velocity along the string for t=t, t=t-dt and t=t-2dt

        self.F_v = np.zeros((len(self.bha.x)+2,3))      # Force along the string for t=t, t=t-dt and t=t-2dt.  Vectorized format
        self.v_v = np.zeros((len(self.bha.x)+2,3))      # Velocity along the string for t=t, t=t-dt and t=t-2dt.  Vectorized format

        self.HS = np.zeros(int(simulation_length/self.dt))              # Hook speed
        self.HS_sp = np.zeros(int(simulation_length/self.dt))           # Hook speed Setpoint
        self.ROP = np.zeros(int(simulation_length/self.dt))             # Rate of Pe
        self.HL = np.zeros(int(simulation_length/self.dt))              # Hook load
        self.WOB = np.zeros(int(simulation_length/self.dt))             # Weight on bit
        
        
        
        self.time = np.arange(0,simulation_length,simulation_length / int(simulation_length/self.dt))

        # Simulation constants
        self.rho = 490.6                            # lb/ft^3
        self.G = 1.671e9                            # lbf/ft^2
        self.E = 2.282e9                            # lbf/ft^2  


        # TD Control System Parameters
        self.kp = 6500                             # p-gain, top drive PI controller
        self.ki = 2500                              # i-gain, top drive PI controller
        self.TD_INERTIA = 750                       # Top drive inertia
        self.ramp = 60/60*2*np.pi                   # maximum set point ramp speed (rpm / second)

        # Initialize controller variables
        self.td_speed_error = 0
        self.td_integrated_speed_error = 0

        # Simulation control
        self.idx = 0                                        # current simulation index
        self.maxsteps = int(simulation_length/self.dt)      # max number of timesteps in simulation
        self.simlength = simulation_length                  # max simulation length in seconds

        # estimation
        self.bit_estimates = np.zeros((2,2))                # allocate storage for bit estimation
        self.bit_constant = bit_constant                            # initial guess
        self.gain = 25                                  # observer gain
        self.update_period = update_period                  # update frequency for ensemble updating




    def step(self, dt=0.1, TD_RPM_SP=None, WOB_SP=None, bit_constant=None):
        """
        Step the simulation forward by dt seconds.  Computations are vectorized in x
        
        Parameters
        ==========
        dt : int
            Number of seconds to advance the simulation
        TD_RPM_SP : float
            A change in the top drive set point. Default is None, meaning no change in top drive set point
        WOB_SP : float
            Current WOB setpoint
        """


        # Check if maximum simulation length has been exceeded.
        # TODO: allow simulation to extend time.
        if self.idx >= self.maxsteps:
            print('Simulation Complete.  Maximum simulation length reached.')
            return 

        # Number of steps to advance the simulation
        nsteps = int(dt/self.dt)

        # To remove high frequency noise from the simulation (to allow for a longer dt and larger dx), 
        # omega nd tau are convolved with the kernel below to remove 40% of energy
        kernel = np.array([0.2, 0.6, 0.2])

        # Loop through time
        for i in range(0,nsteps):

            # advance tau and omega
            self.tau_v[:,:2] = self.tau_v[:,1:]
            self.omega_v[:,:2] = self.omega_v[:,1:]

            # Compute the new setpoint rpm
            if TD_RPM_SP != None:

                # if first iteration of simulation, last rpm setpoint is 0
                if self.idx > 0:
                    last = self.omega_td_sp[self.idx-1]
                else:
                    last = 0

                # increase or decrease the set point at the maximum 'ramp'
                if abs(TD_RPM_SP/60*2*np.pi - last) > self.ramp*self.dt:
                    if last > TD_RPM_SP/60*2*np.pi:
                        self.omega_td_sp[self.idx] = self.omega_td_sp[self.idx-1] - self.ramp*self.dt
                    else:
                        self.omega_td_sp[self.idx] = self.omega_td_sp[self.idx-1] + self.ramp*self.dt
                else:
                    self.omega_td_sp[self.idx] = TD_RPM_SP/60*2*np.pi
           
            # keep the last set point
            else:
                if self.idx > 0:
                    self.omega_td_sp[self.idx] = self.omega_td_sp[self.idx-1]
                else:
                    self.omega_td_sp[self.idx] = 0

            # update the bit constant if provided
            if bit_constant != None:
                self.bit_constant = bit_constant

            #
            # TOP DRIVE
            #

            if self.idx > 0:
                omega_td = self.omega_td[self.idx-1]
                tau_td = self.tau_td[self.idx-1]
            else:
                omega_td = 0
                tau_td = 0

            # Calculate speed error and integrated speed error
            self.td_speed_error = self.omega_td_sp[self.idx] - omega_td
            self.td_integrated_speed_error += self.td_speed_error*self.dt


            # Calculate the new omega and tau of the top drive, assuming a simple PI controller and a motor inertia
            self.tau_td[self.idx] = tau_td + (self.kp * self.td_speed_error + self.ki * self.td_integrated_speed_error) * self.dt
            self.omega_td[self.idx] = omega_td + (self.tau_td[self.idx] - self.tau_v[1,1])/self.TD_INERTIA*self.dt



            # 
            # BIT
            #
            self.tau_bit[self.idx] =  self.bit_constant * self.omega_v[-2,1]
            self.omega_bit[self.idx] = self.omega_v[-2,1]


            # set up the first and last elements of the drillstring arrays for vectorized computation
            self.tau_v[0,1] = self.tau_td[self.idx] - self.TD_INERTIA*(self.omega_td[self.idx] - omega_td) / self.dt
            self.tau_v[-1,1] = self.tau_bit[self.idx]
            self.omega_v[0,1] = self.omega_td[self.idx]
            self.omega_v[-1,1] = self.omega_bit[self.idx]
            

                
            #
            # DRILLSTRING
            #

            # compute the source term for each drillstring element.  The source term consists of:
            #   Inertial damping
            #   Viscous damping
            #   Coulomb friction

            # friction component
            friction = self.Fn*self.bha.OD/12/2

            # inertial damping component
            inertia = self.kt * self.rho*self.bha.J*(self.omega_v[1:-1,1]-self.omega_v[1:-1,0])/self.dt

            # viscous damping component
            viscous = self.ks * self.omega_v[1:-1,1]

            # total source term
            self.source = (friction * self.mu_k + inertia + viscous) / self.dx
            self.source[self.source<0] = -1*self.source[self.source<0]

            # torque at each drillstring element
            self.tau_v[1:-1,2] = self.tau_v[1:-1,1] - self.bha.J*self.G*self.dt/self.bha.l*(self.omega_v[2:,1] - self.omega_v[:-2,1])

            # rpm at each drillstring element
            self.omega_v[1:-1,2] = self.omega_v[1:-1,1] - 2*self.dt/(self.bha.J*self.rho)*((self.tau_v[2:,1] - self.tau_v[:-2,1])/(2*self.bha.l) + self.source)

            # rpm is zero if frictional torque is not exceeded            
            self.omega_v[1:-1,2] = self.omega_v[1:-1,2] * ((np.abs(self.tau_v[1:-1,2]) > friction*self.mu_k) + (np.abs(self.omega_v[1:-1,2]) > self.omega_c)).astype(int)

            self.idx += 1

            # end the simulation at max time
            # TODO: allow extension of simulation time
            if self.idx >= self.maxsteps:
                print('Simulation Complete.  Maximum simulation length reached.')
                return

        # save variables
        self.tau_v[0,2] = self.tau_td[self.idx-1]
        self.tau_v[-1,2] = self.tau_bit[self.idx-1]
        self.omega_v[0,2] = self.omega_td[self.idx-1]
        self.omega_v[-1,2] = self.omega_bit[self.idx-1]

        # Remove 40% of the energy in the system
        self.omega_v[:,1] = np.convolve(self.omega_v[:,1], kernel, mode='same')
        self.omega_v[:,2] = np.convolve(self.omega_v[:,2], kernel, mode='same')
        self.tau_v[:,1] = np.convolve(self.tau_v[:,1], kernel, mode='same')
        self.tau_v[:,2] = np.convolve(self.tau_v[:,2], kernel, mode='same')


    #@jit(nopython=True)
    def speedhelper(sim, dt=0.1, TD_RPM_SP=None, WOB_SP=None, bit_constant=None):
        """
        Step the simulation forward by dt seconds.  Computations are vectorized in x
        
        Parameters
        ==========
        dt : int
            Number of seconds to advance the simulation
        TD_RPM_SP : float
            A change in the top drive set point. Default is None, meaning no change in top drive set point
        WOB_SP : float
            Current WOB setpoint
        """


        # Check if maximum simulation length has been exceeded.
        # TODO: allow simulation to extend time.
        if sim.idx >= sim.maxsteps:
            print('Simulation Complete.  Maximum simulation length reached.')
            return 

        # Number of steps to advance the simulation
        nsteps = int(dt/sim.dt)

        # To remove high frequency noise from the simulation (to allow for a longer dt and larger dx), 
        # omega nd tau are convolved with the kernel below to remove 40% of energy
        kernel = np.array([0.2, 0.6, 0.2])

        # Loop through time
        for i in range(0,nsteps):

            # advance tau and omega
            sim.tau_v[:,:2] = sim.tau_v[:,1:]
            sim.omega_v[:,:2] = sim.omega_v[:,1:]

            # Compute the new setpoint rpm
            if TD_RPM_SP != None:

                # if first iteration of simulation, last rpm setpoint is 0
                if sim.idx > 0:
                    last = sim.omega_td_sp[sim.idx-1]
                else:
                    last = 0

                # increase or decrease the set point at the maximum 'ramp'
                if abs(TD_RPM_SP/60*2*np.pi - last) > sim.ramp*sim.dt:
                    if last > TD_RPM_SP/60*2*np.pi:
                        sim.omega_td_sp[sim.idx] = sim.omega_td_sp[sim.idx-1] - sim.ramp*sim.dt
                    else:
                        sim.omega_td_sp[sim.idx] = sim.omega_td_sp[sim.idx-1] + sim.ramp*sim.dt
                else:
                    sim.omega_td_sp[sim.idx] = TD_RPM_SP/60*2*np.pi
           
            # keep the last set point
            else:
                if sim.idx > 0:
                    sim.omega_td_sp[sim.idx] = sim.omega_td_sp[sim.idx-1]
                else:
                    sim.omega_td_sp[sim.idx] = 0

            # update the bit constant if provided
            if bit_constant != None:
                sim.bit_constant = bit_constant

            #
            # TOP DRIVE
            #

            if sim.idx > 0:
                omega_td = sim.omega_td[sim.idx-1]
                tau_td = sim.tau_td[sim.idx-1]
            else:
                omega_td = 0
                tau_td = 0

            # Calculate speed error and integrated speed error
            sim.td_speed_error = sim.omega_td_sp[sim.idx] - omega_td
            sim.td_integrated_speed_error += sim.td_speed_error*sim.dt


            # Calculate the new omega and tau of the top drive, assuming a simple PI controller and a motor inertia
            sim.tau_td[sim.idx] = tau_td + (sim.kp * sim.td_speed_error + sim.ki * sim.td_integrated_speed_error) * sim.dt
            sim.omega_td[sim.idx] = omega_td + (sim.tau_td[sim.idx] - sim.tau_v[1,1])/sim.TD_INERTIA*sim.dt



            # 
            # BIT
            #
            sim.tau_bit[sim.idx] =  sim.bit_constant * sim.omega_v[-2,1]
            sim.omega_bit[sim.idx] = sim.omega_v[-2,1]


            # set up the first and last elements of the drillstring arrays for vectorized computation
            sim.tau_v[0,1] = sim.tau_td[sim.idx] - sim.TD_INERTIA*(sim.omega_td[sim.idx] - omega_td) / sim.dt
            sim.tau_v[-1,1] = sim.tau_bit[sim.idx]
            sim.omega_v[0,1] = sim.omega_td[sim.idx]
            sim.omega_v[-1,1] = sim.omega_bit[sim.idx]
            

                
            #
            # DRILLSTRING
            #

            # compute the source term for each drillstring element.  The source term consists of:
            #   Inertial damping
            #   Viscous damping
            #   Coulomb friction

            # friction component
            friction = sim.Fn*sim.bha.OD/12/2

            # inertial damping component
            inertia = sim.kt * sim.rho*sim.bha.J*(sim.omega_v[1:-1,1]-sim.omega_v[1:-1,0])/sim.dt

            # viscous damping component
            viscous = sim.ks * sim.omega_v[1:-1,1]

            # total source term
            sim.source = (friction * sim.mu_k + inertia + viscous) / sim.dx
            sim.source[sim.source<0] = -1*sim.source[sim.source<0]

            # torque at each drillstring element
            sim.tau_v[1:-1,2] = sim.tau_v[1:-1,1] - sim.bha.J*sim.G*sim.dt/sim.bha.l*(sim.omega_v[2:,1] - sim.omega_v[:-2,1])

            # rpm at each drillstring element
            sim.omega_v[1:-1,2] = sim.omega_v[1:-1,1] - 2*sim.dt/(sim.bha.J*sim.rho)*((sim.tau_v[2:,1] - sim.tau_v[:-2,1])/(2*sim.bha.l) + sim.source)

            # rpm is zero if frictional torque is not exceeded            
            sim.omega_v[1:-1,2] = sim.omega_v[1:-1,2] * ((np.abs(sim.tau_v[1:-1,2]) > friction*sim.mu_k) + (np.abs(sim.omega_v[1:-1,2]) > sim.omega_c)).astype(int)

            sim.idx += 1

            # end the simulation at max time
            # TODO: allow extension of simulation time
            if sim.idx >= sim.maxsteps:
                print('Simulation Complete.  Maximum simulation length reached.')
                return

        # save variables
        sim.tau_v[0,2] = sim.tau_td[sim.idx-1]
        sim.tau_v[-1,2] = sim.tau_bit[sim.idx-1]
        sim.omega_v[0,2] = sim.omega_td[sim.idx-1]
        sim.omega_v[-1,2] = sim.omega_bit[sim.idx-1]

        # Remove 40% of the energy in the system
        sim.omega_v[:,1] = np.convolve(sim.omega_v[:,1], kernel, mode='same')
        sim.omega_v[:,2] = np.convolve(sim.omega_v[:,2], kernel, mode='same')
        sim.tau_v[:,1] = np.convolve(sim.tau_v[:,1], kernel, mode='same')
        sim.tau_v[:,2] = np.convolve(sim.tau_v[:,2], kernel, mode='same')

        return sim


    def initialize_estimator(self, array_size, intensity):
        """
        Initialize the estimator.  
        
        Parameters
        ==========
        array_size : int
            Length to preallocate for estimation

        intensity : float
            Intensity of update
        
        """
        self.bit_estimates = np.zeros((array_size,2))
        self.bit_estimates[:] = np.nan
        self.estimator_idx = 0
        self.error = 0
        self.intensity = intensity
        self.bit_constant0 = self.bit_constant
        self.t0 = 0

    def update_estimate(self, tau_act):
        """
        Update the parameter estimates based on the actual value of surface torque

        Parameters
        ==========
        tau_act : Float
            the actual surface torque
        """

        if self.estimator_idx >= len(self.bit_estimates):
            return
        err = tau_act - self.tau_td[self.idx-1]
        self.bit_constant += self.gain * err * self.dt # * (np.random.random_sample())

        
        #self.bit_constant = self.bit_constant0 + self.intensity * expit((self.time[self.idx-1] - self.t0 - 6)/ 12)
        self.bit_constant = max(0, self.bit_constant)
        self.bit_estimates[self.estimator_idx,:] = [self.time[self.idx-1], self.bit_constant]
        self.estimator_idx += 1
        self.error += np.sqrt(tau_act ** 2 - self.tau_td[self.idx-1] ** 2)


    def clone_estimator(self, source):
        """
        For ensemble based estimation, clone self to match best estimator

        Parameters
        ==========
        source : DrillPyze.Simulation
            the best estimator
        """

        try:

            self.bit_constant = source.bit_constant
            self.bit_constant0 = source.bit_constant
            self.t0 = self.time[self.idx-1]
            self.tau_v = np.copy(source.tau_v)
            self.omega_v = np.copy(source.omega_v)

            self.tau_td[self.idx] = source.tau_td[self.idx]
            self.tau_bit[self.idx] = source.tau_bit[self.idx]
            self.omega_td[self.idx] = source.omega_td[self.idx]
            self.omega_bit[self.idx] = source.omega_bit[self.idx]

        except:
            return -1
        
        self.error = 0

    def __lt__(self, comp):
        """
        Redefining the less than method to allow for sorting.  One simulation is less than
        another simulation iff estimation error over a specified timeframe is smaller.  The
        timeframe is defined when the simulations are initialized.

        This method does not check for consistency.

        Parameters
        ==========
        comp : DrillPyze.Simulation
            the simulation to compare against
        """

        return self.error < comp.error
    
    def parstep(args):
        """
        Helper function for multiprocessing
        """
        Simulation.speedhelper(args[0], args[1][0], args[1][1])
    


class Estimator:
    """
    A class to maintain and update an ensemble based estimation technique
    """



    def __init__(self, num_instances, wellbore, bha):
        """
        Initialize an pool of workers (only a ThreadingPool works currently)
        to create a series of potential estimates

        Parameters
        ==========
        num_instances : int
            Number of instances or workers for the ensemble.  Assumes that the
            number of workers and the number of estimates is equal
        wellbore : DrillPyze.well.wellbore
            Wellbore for the simulation
        bha : DrillPyze.drillstring.BHA
            BHA to be used for the simulation
        """


        # Create an array of simulations
        self.simulations = [None] * num_instances
        self.num_sim = num_instances
        for i in range(num_instances):
            self.simulations[i] = Simulation(i, wellbore, bha)

        # Create a multiprocessing pool
        self.pool = Pool(nodes=num_instances)           # set up the multiprocessing pool to increment the simulations


    def __del__(self):
        """
        Ensure the processing pool is closed when the Estimator exits
        """
        self.pool.terminate()
        self.pool.join()


    def initialize(self, bit_depth, sim_length, estimate_freq):
        """
        Initialize the simulations to the same start point

        Parameters
        ==========
        bit_depth : float
            Current bit depth
        sim_length : int
            Number of seconds to simulate
        estimate_freq : float
            Frequency (in Hz) of available data
        """
        for i in range(self.num_sim):
            # self.simulations[i].initialize_simulation(bit_depth, sim_length, bit_constant=((i+1)/3*1000))
            self.simulations[i].initialize_simulation(bit_depth, sim_length)
            self.simulations[i].initialize_estimator(int(sim_length/estimate_freq), ((i-self.num_sim/2)/self.num_sim)*500)



    def increment(self, time, TD_SETPOINT=0, TD_ACT=0):
        """
        Increment the simulations with a certain topdrive rpm and recorded top
        drive torque.  The simulations simulate forward at the topdrive rpm
        and then update their estimation assuming the topdrive torque is measured
        at the end of the simulation increment.

        The simulation pool is then sorted by smallest error and all simulations
        are cloned to be equal to the best simulation.

        Parameters
        ==========
        time : float
            Number of seconds to step the simulations forward
        TD_SETPOINT : float
            Top drive rpm setpoint.  Default is 0
        TD_ACT : float
            Top drive actual torque.  Default is 0
        """

        # Activate the pool
        self.pool.restart()

        # Create the argument tuples for each simulation
        args = [(time, TD_SETPOINT)] * self.num_sim

        # Map the simulations and arguments to the workers and start the pool
        self.pool.map(Simulation.parstep, zip(self.simulations, args))
        self.pool.close()

        # wait for the pool to finish
        self.pool.join()

        # update the estimates for each worker
        for i in range(self.num_sim):
            self.simulations[i].update_estimate(TD_ACT)

        # Sort the simulations based on lowest error
        self.simulations.sort()

        # Clone the best estimate
        for i in range(1,self.num_sim):
            self.simulations[i].clone_estimator(self.simulations[0])



