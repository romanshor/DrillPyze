if __name__ == '__main__':

    from DrillPyze import well, drillstring, statics, dynamics
    import numpy as np
    import matplotlib.pyplot as plt
    import time

    from pathos.helpers import freeze_support, shutdown
    freeze_support()

    sim_length = 60

    well1 = well.wellbore({'xlsx':'../data/survey.xlsx'}, 1000, 5)
    bha1 = drillstring.BHA('BHA 1', '../data/bha.xlsx')



    plot_frequency = 0.5

    print("Initializing Estimation Pool")
    estimator = dynamics.Estimator(10, well1, bha1)
    estimator.initialize(7000, sim_length, plot_frequency)





    print("Iterating")
    tic = time.perf_counter()
    for i in range (0,int(sim_length*plot_frequency)):
        if i > 60*plot_frequency:
            estimator.increment(1/plot_frequency, 60, 21000)
        elif i > 1*plot_frequency:
            estimator.increment(1/plot_frequency, 30, 24000)
        else:
            estimator.increment(1/plot_frequency)

    toc = time.perf_counter()
    print(f"Estimation completed in {toc - tic:0.4f} seconds") 


    plt.ion()
    fig = plt.figure()
    ax1 = fig.add_subplot(3,4,(3,4))
    ax2 = fig.add_subplot(3,4,(7,8))
    ax3 = fig.add_subplot(3,4,(1,5))
    ax4 = fig.add_subplot(3,4,(2,6))
    ax5 = fig.add_subplot(3,4,(11,12))

    

    for sim in estimator.simulations:
        

        ax1.plot(sim.time, sim.omega_bit)
        ax1.plot(sim.time, sim.omega_td, 'b')
        ax1.plot(sim.time, sim.omega_td_sp, 'k')


        ax2.plot(sim.time, sim.tau_td, 'b')
        ax2.plot(sim.time, sim.tau_bit, 'r')


        ax3.plot(sim.tau_v[1:-1,2],bha1.x,'b')

        #ax4.plot(sim.omega[:,2]/(np.pi)*60,bha1.x, 'b')
        ax4.plot(sim.omega_v[1:-1,1],bha1.x,'b')

        ax5.plot(sim.bit_estimates[:,0], sim.bit_estimates[:,1])

        fig.canvas.draw()
        fig.show()


    shutdown()

if __name__ == '__del__':
    from pathos.helpers import shutdown
    shutdown()


