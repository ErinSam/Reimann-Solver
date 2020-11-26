"""
    Python script that implements Total Variation Diminishing (TVD) on the MUSCL-Hancock
    scheme to solve 1D Euler Equations
"""
# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electrical & Electronics Engg. Minor



import numpy as np 
import math as m
import matplotlib.pyplot as plt
import slope_limiters as slp
import compressed_erin as cmprin






def intialisation(U):
    """
        Function that initialises the initial data line

        Args:
            U: ndarray(M+2,3); the conserved elements 2D array for the time step 

        Returns:
            U: ndarray(M+2,3); the conserved elements 2D array for the time step after applying 
                                the initial coniditons to it
    """
    # TODO



def boundary_conditions(U):
    """
        Function that applies the boundary conditions according to our choice of 
        boundary conditions
        
        Args: 
            U: ndarray(M+2,3); the conserved elements 2D array for the time step 

        Returns: 
            U: ndarray(M+2,3); the conserved elements 2D array for the time step after applying 
                                the boundary coniditons to it
    """
    # Transmissive Boundary Conditions
    U[0], U[M+1] = U[1], U[M]

    # Relfective Boundary Conditions 
    # TODO

    return U



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



def time_step_size(U, cfl, dx):
    """ 
        Function that calculates the size of the time step

        Args:
            U: ndarray(M+2,3); the conserved elements 2D array for the time step 
            cfl: float; courant number 
            dx: float; mesh size

        Returns:
            dt: float; the size of the time step 
    """
    c_ratio = 1.4
    
    conv_prim = np.vectorize(cmprin.consv_prim_1D)
    W = conv_prim(U)

    speeds = abs(W[:,1]) + m.sqrt(c_ratio * W[:,2] / W[:,0])

    dt = cfl * dx / np.max(speeds)

    return dt



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

    # Obtain the following from the user
    print("\n\nEnter the following data")
    cfl = float(input("CFL Number: "))
    num_data_pts = int(input("Number of data points: "))
    time = float(input("Time to which the solution should be obtained: "))
    print("\nThe length of the data line is assumed to be 1m. Starts from 0 and ends at 1m")

    # Creating the arrays given the user input 
    U = np.zeros((M+2, 3))

    # Initialising the array
    # TODO

    # Applying the boundary conditions
    U = boundary_conditions(U)

    # Calling the iteration step (marching in time)
    # TODO







if __name__ == "__main__":
    main()

