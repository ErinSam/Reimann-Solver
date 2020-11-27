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



def bev(U, dt, dx, M):
    """ BOUNDARY EXTRAPOLATED VALUES
        Function that calculates the boundary extrapolated values. Results are used to obtain the
        evolved BEVs and they depend a lot on the choice of the slope limter chosen.

        Args:
            U: ndarray(M+2,3); the conserved elements 2D array for the time step 
            dt: float; size of the time step 
            dx: float; mesh size
            M: int; number of mesh points

        Returns:
    """
    # Creating ndarray for boundary extrapolated values
    U_bev = np.zeros((M+2,2,3))

    # Obtaining the intercell element-wise jump vector (intrcll_delta)
    intrcll_delta = U[1:] - U[:-1]

    # Calculating upwind ratio
    r = intrcll_delta[:-1] / intrcll_delta[1:]

    # Calculating the initial slopes
    delta = 0.5 * (1+omega) * intrcll_delta[:-1] + 0.5 * (1-omega) * intrcll_delta[1:]

    # Obtaining slope limiter
    XI = np.zeros(r.shape)
    for i, val in enumerate(r):
        XI[i] = slp.SUPERBEE(val)                       # FIXME

    # Calculcating the limited slopes
    limited_slope = XI * delta
    
    # Obtaining boundary extrapolated values
    U_bev[1:-1,0] = U[1:-1] - limited_slope/2
    U_bev[1:-1,1] = U[1:-1] + limited_slope/2

    # Applying zero slope boundary conditions to the boundary BEVs
    U_bev[0,:] = U[0]
    U_bev[-1,:] = U[-1]

    # Obtaining Evolved BEVs
    U_ebev = evolved_bev(U_bev, dt, dx)

    return U_ebev



def evolved_bev(U_bev, dt, dx):
    """
        Function that caluclates the evolved boundary extrapolated values 
        This results of this function are then used to obtain the intercell flux by solvng the RP

        Args:
            U_bev: ndarray(M+2,2,3); the boundary extrapolated values of the conserved elements
            dt: float; size of the time step 
            dx: float; mesh size

        Returns:
            U_ebev: ndarray(M+2,2,3); evolved boundary extrapolated values of the conserved vars
    """
    # Creatin ndarray for eveolved boundary extrapolated values
    U_ebev = np.zeros(U_bev.shape)

    # Calculating the flux function
    F = np.empty(U_bev.shape)
    for i, val_1 in enumerate(U_bev):
        for j, val_2 in enumerate(val_1):
            F[i,j] = cmprin.flux_euler_1D(val_2)

    # Calculating evolved boundary extrapolated values
    U_ebev[:,0,:] = U_bev[:,0,:] + 0.5*dt/dx * (F[:,0,:] - F[:,1,:])
    U_ebev[:,1,:] = U_bev[:,1,:] + 0.5*dt/dx * (F[:,0,:] - F[:,1,:])

    return U_ebev 



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

