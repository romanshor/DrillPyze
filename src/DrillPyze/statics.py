import pandas as pd
import numpy as np
from DrillPyze import *


def normal_forces(well, bha, WOB=0, TOB=0):
    """
    Calculate the normal forces along the drillstring.

    Reference:
    Johancsik, Friesen and Dawson. Torque and Drag in Directional Wells – Prediction and Measurement.  Journal of Petroleum Technology.  1984


    Parameters
    ==========
    well : DrillPyze.well.well_bore
        A well_bore object initialized with a well profile.  The wellpath must already be interpolated.
    bha : DrillPyze.drillstring.BHA
        A BHA object initialized with a drillstring configuration.  The BHA must already be discretized.
    WOB : float
        axial force at the bit (weight on bit).  Default value is 0
    TOB : float
        torque at the bit.  Default value = 0
    """

    n_elem = len(bha.x)
    mu = 0.3

    Ft = np.zeros(n_elem)       # Tension
    Ft_pu = np.zeros(n_elem)    # Tension, pulling up
    Ft_so = np.zeros(n_elem)    # Tension, slacking off
    Fn = np.zeros(n_elem)       # Normal force
    M = np.zeros(n_elem)        # Torque


    # convert to radians
    azi = well.azi / 180 * np.pi
    inc = well.inc / 180 * np.pi

    # Calculate the change in azimuth and inclination across each elements
    d_azi = azi[1:] - azi[:-1]
    d_inc = inc[1:] - inc[:-1]

   

    for i in range(0,n_elem):

        if i==0:
            # Boundary condition at the bit
            M[n_elem-1] = TOB
            Ft[n_elem-1] = WOB
        else:
            # Calculate normal force for the element.  Eqn 1
            fn = np.sqrt(np.power(Ft[n_elem-i]*d_azi[n_elem-1-i]*np.sin(inc[n_elem-i]), 2) +
                np.power(Ft[n_elem-i]*d_inc[n_elem-1-i] + bha.w[n_elem-i]*np.sin(inc[n_elem-i]), 2)
                )
            # Calculate the change in tension.  Eqn 2
            dft = bha.w[n_elem-i]*np.cos(inc[n_elem-i])
            # Calcualte the friction formce
            fric = mu * fn

            # Calculate the change in torque.  Eqn 3
            dM = mu * fn * bha.OD[n_elem-i]/2/12

            # Calculate and store the tension, torque and normal forces
            Ft[n_elem-1-i] = Ft[n_elem-i] + dft
            Ft_pu[n_elem-1-i] = Ft_pu[n_elem-i] + dft + fric
            Ft_so[n_elem-1-i] = Ft_so[n_elem-i] + dft - fric
            Fn[n_elem-i] = fn
            M[n_elem-1-i] = M[n_elem-i] + dM
    
    # kernel = np.array([0.1, 0.8, 0.1])
    # Fn[1:-1] = np.convolve(Fn, kernel, mode='same')[1:-1]
    # M[1:-1] = np.convolve(M, kernel, mode='same')[1:-1]


    return (Ft, Ft_pu, Ft_so, Fn, M)





