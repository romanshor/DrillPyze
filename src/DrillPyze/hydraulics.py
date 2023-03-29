"""
Hydraulics module for DrillPyze.

This module contains functions for calculating the flow rate, fluid velocity,
Reynolds number for a pipe and the frictional pressure for a pipe and an annulus.
"""

def fluid_velocity(Q, Do, Di=None):
    """
    Calculate the fluid velocity in a pipe or an annulus.
    
    Parameters
    ==========
    Q : float
        Flow rate in gpm
    Do : float
        Outer diameter of the pipe in inches
    Di : float
        Inner diameter of the pipe in inches.
        If Di is None, then the function will calculate the fluid velocity for a pipe.
    """

    # For pipe flow
    if Di is None:
        v = Q / (2.448 * Do**2)
    
    # For annular flow 
    else:
        v = Q / (2.448 * (Do**2 - Di**2))

    # Return the fluid velocity
    return v


def reynolds_number(rho, v, mu, Do, Di=None):
    """
    Calculate the Reynolds number for a pipe or an annulus.

    Parameters
    ==========
    rho : float
        Fluid density in lbm/ft^3
    v : float
        Fluid velocity in ft/s
    mu : float
        Fluid viscosity in cp
    Do : float
        Outer diameter of the pipe in inches
    Di : float
        Inner diameter of the pipe in inches.
        If Di is None, then the function will calculate the Reynolds number for a pipe.
    """

    # For pipe flow
    if Di is None:
        NRe = 928 * rho * v * Do / mu
    
    # For annular flow
    else:
        NRe = 757 * rho * v * (Do - Di) / mu

    # Return the Reynolds number
    return NRe

def friction_factor(NRe):
    """
    Calculate the friction factor for a pipe or an annulus.

    Parameters
    ==========
    NRe : float
        Reynolds number for the pipe
    """

    # If in laminar flow, use the laminar flow equation
    if NRe < 2000:
        f = 16 / NRe

    # If in turbulent flow, use the turbulent flow equation
    else:
        f = 0.0791 / (NRe ** (1/4))

    # Return the friction factor
    return f


def frictional_pressure_drop(rho, Q, mu, L, Do, Di=None):
    """
    Calculate the frictional pressure drop in a pipe or an annulus.

    Parameters
    ==========
    rho : float
        Fluid density in ppg
    Q : float
        Flow rate in gpm
    mu : float
        Fluid viscosity in cp
    L : float
        Length of the pipe in ft
    Do : float
        Outer diameter of the pipe in inches
    Di : float
        Inner diameter of the pipe in inches.
        If Di is None, then the frictional pressure drop is calculated for a pipe.
    """

    # Calculate fluid velocity
    v = fluid_velocity(Q, Do, Di)

    # Calculate Reynolds number
    NRe = reynolds_number(rho, v, mu, Do, Di)

    # If in laminar flow, use the laminar flow equation
    if NRe < 2000:
        if Di is None:
            dP = mu * v / (1500 * Do**2) * L
        else:
            dP = mu * v / (1000 * (Do**2 - Di**2)) * L

    # If in turbulent flow, use the turbulent flow equation
    else:
        if Di is None:
            dP = rho**0.75*v**1.75*mu**0.25/(1800 * Do**1.25)*L
        else:
            dP = rho ** 0.75 * v ** 1.75 * mu ** 0.25 / (1800 * (Do - Di) ** 1.25) * L
    
    # Return the pressure drop
    return dP