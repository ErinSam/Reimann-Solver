"""
    Python scripts that applies Godunov Method to solve a system of linear, hyperbolic
    equations

    System of Equations is of the form:
        W_t + F(W)_x = 0
            
            where,
                W_t is the partial derivatives of W wrt t 
                W_x is the partial derivative of W wrt x 
                F is the flux function
    
    Boundary Condition: Transmissive Boundary Conditons
"""
# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electrical & Electronics Engg. Minor 



import numpy as np 
import math as m 
import scipy.integrate as integrate
from iterative_scheme_for_pressure import *



def initialisation(W, dx):
    """ 
        Function that initialises the Primitive Variable vector W at t = 0 according 
        to their definitions using the finite volume method. 
        
        Args:
            W - ndarray(3, num_data_pts + 2) ; the primitive variable vector 

        Returns:
            W - ndarray(3, num_data_pts + 2): the initialised variable vector
    """
    
    def rho_init():
        return rho
    
    def u_init():
        return u 

    def p_init():
        return p

    # Creating a list of initialisation functions
    # each element in the list is the initiaisation function of a primitive variable
    init_func_list = [rho_init, u_init, p_init]

    # Initialising the values at t = 0
    # Note that we are using the finite volume method 
    for init, i in enumerate(init_func_list):
        for prim_var, j in enumerate(W[i,1:-1]):
           prim_var = (1/dx) * integrate.quad(init, j-0.5, j+0.5)
    
    W[:,0], W[:,-1] = W[:,1], W[:,-2]

    return W
    

def reimann_problem(W_L, W_R):
    """
        Function that solves the Reimann Problem to obtain the values at the intercell 
        boundary and also provide largest wave speed 
        
        Args:
            W_L - ndarray(3,) ; the values of the left region in Reimann Problem 
            W_R - ndarray(3,) ; the values of the right region in Reimann Problem 

        Returns: 
            W_0 - ndarray(3,) ; solution of the Reimann Problem at the point x/t = 0
            S_RP_max - float ; absolute value of the largest wave speed that is generated due 
                to the conditions of the RP 
    """

    # INITIALISING VARIABLES
    # Specific Heat Ratio 
    c_ratio = 1.4

    # Denoting the driver region with L and the driven region with R 
    # Obtaining initial conditons for pressure, density and velocity of the driver 
    # and driven region
    rho_L, u_L, p_L = tuple(W_L)
    rho_R, u_R, p_R = tuple(W_R)


    # Speed of sound in the driver and driven regions
    a_L = m.sqrt(c_ratio * p_L/rho_L)
    a_R = m.sqrt(c_ratio * p_R/rho_R)

    # Initialising the values of the data-dependent constants A_L, B_L, A_R, B_R
    A_L = 2 / ((c_ratio+1)*rho_L)
    B_L = p_L * (c_ratio-1)/(c_ratio+1)
    A_R = 2 / ((c_ratio+1)*rho_R)
    B_R = p_R * (c_ratio-1)/(c_ratio+1)


    # Obtaining the Pressure and Velocity of the Star Region 
    p_star, u_star = star_region_properties(p_L=p_L, p_R=p_R, 
                                            rho_L=rho_L, rho_R=rho_R, 
                                            u_L=u_L, u_R=u_R, 
                                            a_L=a_L, a_R=a_R, 
                                            A_L=A_L, A_R=A_R, 
                                            B_L=B_L, B_R=B_R, 
                                            c_ratio=c_ratio)   





    return W_0, S_RP_max


def godunov_method_sys_nonlinear_hyperbolic():
    """
        Function that uses Godunov's Method's second version to obtain a solution of the 
        system of nonlinear, hyperbolic equations over the temporal space
        
        User defines the length of the spatial domain of interest and the number of data
        points to sample data for 
        
        Plots are produced at the end 
    """

    # Obtaining user's input
    length = float(input("\nEnter the length of the data line: ")) 
    num_data_pts = int(input("Enter the number of data (the precision to which) the 
                                solution should be obtained for: "))
    time_period = float(input("\nEnter the time period for which the solution should
                                be obtained"))
    cfl = float(input("CFL coefficient: "))
    print("\nUser should note that the value incremental value for each time step will 
                                be obtained by applying the CFL stability condition."))

    # Creating Primitive Variable Vector, W
    W = np.zeros((3, num_data_pts + 2))
    
    # Initialising the values of W 
    W = initialisation(W)
 







 
