"""
    Python script that implements Total Variation Diminishing (TVD) on the MUSCL-Hancock
    scheme to solve 1D Euler Equations
"""
# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electrical & Electronics Engg. Minor



import numpy as np 
import math as m
import matplotlib.pyplot as plt



if __name__ == "__main__":
    main()



def intialisation():
    """
        Function that initialises the initial data line

        Args:

        Returns:
    """
    # TODO



def boundary_conditions():
    """
        Function that applies the boundary conditions according to our choice of 
        boundary conditions
        
        Args: 

        Returns: 
    """
    # Transmissive Boundary Conditions
    # TODO

    # Relfective Boundary Conditions 
    # TODO



def bev():
    """ BOUNDARY EXTRAPOLATED VALUES
        Function that calculates the boundary extrapolated values. Results are used to obtain the
        evolved BEVs and they depend a lot on the choice of the slope limter chosen.

        Args:

        Returns:
    """
    # TODO



def evolved_bev():
    """
        Function that caluclates the evolved boundary extrapolated values 
        This results of this function are then used to obtain the intercell flux by solvng the RP

        Args:

        Returns:
    """
    # TODO


def RP_solver():
    """
        Function that calculates value of the intercell numerical flux at the origin of the RP
        The RP(bar_U_L, bar_U_R)

        Args:

        Returns:
    """
    # TODO



def time_step_size():
    """ 
        Function that calculates the size of the time step

        Args:

        Returns:
    """
    # TODO



def stepping_time():
    """ 
        Function that calculates the values of the conserved variables at the next time step

        Args:

        Returns:
    """
    # TODO



def plotter()
    """ 
        Plots stuff

        Args:

        Returns:
    """
    # TODO



def main():
    # TODO




