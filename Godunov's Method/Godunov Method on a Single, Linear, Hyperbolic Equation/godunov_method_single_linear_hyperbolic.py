"""
    Python scripts that applies the Godunov Method to solve a Single, Linear, 
    Hyperbolic Equation
"""
# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electronics & Electrical Engg. Minor


import numpy as np 
import scipy.integrate as integrate
import math as m 



def reimann_problem(rp):
    """
        Function that provides an exact solution to a Reimann Problem at x/t = 0

        Args:
            rp - ndarray(2,) ; list that contains the initial values of the reimann 
                problem for the left and right sides 

        Returns:
            u_half - float ; exact solution to the Reimann Problem defined by u at 
                x/t = 0
            S - float ; maximum characteristic speed in Reimann Problem
                        in case shock is produced, it is equal to the shock speed
                        in case rarefaction is produced, 
                            it is equal maximum absolute value of head & tail speed
    """

    if ( rp[0] > rp[1] ):
        shock_speed = 0.5 * ( rp[0] + rp[1] )
        if ( shock_speed >= 0 ):
            return rp[0], abs(shock_speed)
        else: 
            return rp[1], abs(shock_speed)
    else:
        if ( rp[0] >= 0 ):
            return rp[0], max(abs(rp[0]), abs(rp[1]))
        elif ( (rp[0] < 0) & (rp[1] > 0) ):
            return 0, max(abs(rp[0]), abs(rp[1]))
       else: 
            return rp[1], max(abs(rp[0]), abs(rp[1])) 


def flux_function(x):
    """
        Function that is the definition of the flux function our single, linear,
        hyperbolic equation 

        Args: 
            x - float ; value at which we want to evaluate the flux function 

        Returns: 
            f_u - float ; value of the flux function at x
    """
    
    # User can change the definition of the flux function to fit their problem at hand 
    f_u = # ADD definition of flux function 
    
    return f_u


def initiali_tsar(length, num_data_pts):
    """
        Function that initialises the values for Godunov's Method at t = 0

        Args:
            length - float ; length of the data line 
            num_data_pts - int ; number of data points we want to sample for

        Return: 
            u - ndarray(num_data_pts + 2,) ; contains the average value 
    """ 

    def init_func(x):
        val =                                                                                                                           # ADD definition of IVP 
        return val
    
    dx = length/(num_data_pts + 1)
    
    # Creating array to store values of u     
    u = np.zeros(num_data_pts + 2)
    
    # Initialising the values at t = 0
    # Note that we are using the finite volume method
    for j in range(1, num_data_pts + 1):
        u[j] = (1/dx) * integrate.quad(init_func, j-0.5, j+0.5)
    
    u[0], u[-1] = u[1], u[-2]

    return u



def godunov_method_single_linear_hyperbolic():    
    """ 
        This is our main function that uses Godunov's Method to obtain a solution of 
        the single, linear, hyperbolic equation over the time interval specified by the 
        user. 
        The user also provides the length of the region to obtain solutions for and 
        the number of data points (the precision to which) he wants the solution to be 
        obtained for
        This function also produces plots after a defined number of time steps
    """

    length = float(input("\nEnter the length of the data line: ")                                                                       # CHANGE print statement 
    num_data_pts = int(input("Enter the number of data (the precision to which) the 
                                solution should be obtained for: ")
    time_period = float(input("\nEnter the time period for which the solution should
                                be obtained")
    print("\nUser should note that the value incremental value for each time step will 
                                be obtained by applying the CFL stability condition.")
    

    # Initialising the values for Godunov's Method at t = 0
    u = initiali_tsar(length, num_data_pts)

    # Obtaining the value of dt using the CFL stability condition
    dt =                                                                                                                                # ADD function to calc dt  
    num_time_steps = time_period / dt + 1

    # GODUNOV'S METHOD
    u_half = np.zeros(num_data_pts + 1)
    S = np.zeros(num_data_pts + 1)
    flux = np.zeros(num_data_pts + 1)

    for t in range (num_time_steps):
        for i range( num_time_steps + 1):
            rp = np.array([u[i], u[i+1]])
            u_half[i], S[i] = reimann_problem(rp)
            flux[i] = flux_function(u_half[i])
            
        # The Godunov Update Step
        u[1:-1] += np.max(S) / dx * (flux[:-1] - flux[1:])
        
        # Applying the boundary condition 
        u[0], u[-1] = u[1], u[-2]
        
        # ADD plotting feature

