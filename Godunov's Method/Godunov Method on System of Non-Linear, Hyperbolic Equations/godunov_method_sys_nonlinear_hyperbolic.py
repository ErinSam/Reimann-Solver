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
import matplotlib.pyplot as plt



def initialisation(W, dx):
    """ 
        Function that initialises the Primitive Variable vector W at t = 0 according 
        to their definitions using the finite volume method. 
        
        Args:
            W - ndarray(3, num_data_pts + 2) ; the primitive variable vector 

        Returns:
            W - ndarray(3, num_data_pts + 2): the initialised variable vector
    """
    
    def rho_init(x):
        """User must edit this function according to initilisation of density 
            Args:
                x - float; 
            Returns: 
                rho - float; density at (x,0)
        """
        # ADD definition 
        return rho
    
    def u_init(x):
        """User must edit this function according to initialisation of velocity
            Args:
                x - float;
            Returns:
                u - float; velocity at (x,0)
        """
        # ADD definition
        return u 

    def p_init(x):
        """User must edit this function according to initialisation of pressure
            Args:
                x - float;
            Returns: 
                p - float; pressure at (x,0)
        """
        # ADD definition
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


    # OBTAINING WAVE SPEEDS & SOLUTION OF W @ x/t = 0 (W_0)
    # Creating wave_speeds, ndarray, that stores absolute values of wave speeds
    wave_speeds = np.zeros(0)
    # Initialsing W_0 with W_R in case of the trivial case where W_L = W_R
    W_0 = W_R 

    # Obtaining Exact Solution for Left of Contact Discontinuity
    if ( p_star <= p_L ):
        # Left wave is an expansion wave 
        rho_star_L = rho_L * pow(p_star/p_L, 1/c_ratio)
        S_HL = u_L - a_L
        S_TL = u_star - m.sqrt(c_ratio * p_star / rho_star_L)

        # Adding absolute values of wave speeds to wave_speeds
        np.append(wave_speeds, np.abs([S_HL, S_TL])
        
        # Setting W_0 (if possible)
        if ( u_star > 0 ):
            if ( S_TL > 0 ):
                W_0 = W_L
            elif ( S_HL > 0 ):
                W_0[0], W_0[1], W_0[2] = tuple(rho_star_L, u_star, p_star)
            else:
                W_0[0] = p_L * pow(( (2/(c_ratio+1)) 
                            + (c_ratio-1)/(a_L*(c_ratio+1)) * (u_L) ), 2*c_ratio/(c_ratio-1) )
                W_0[1] = 2/(c_ratio+1) * ( a_L + u_L * (c_ratio-1)/2 )
                W_0[2] = rho_L * pow(( (2/(c_ratio+1)) 
                            + (c_ratio-1)/(a_L*(c_ratio+1)) * (u_L) ), 2/(c_ratio-1) )
    else:
        # Left wave is a shock wave
        rho_star_L = rho_L * ( (p_star/p_L) 
                        + (c_ratio-1)/(c_ratio+1) ) / ( (c_ratio-1)/(c_ratio+1)*(p_star/p_L) + 1 )
        S_L = u_L - a_L * m.sqrt( (c_ratio+1)/(2*c_ratio) * p_star/p_L + (c_ratio-1)/(2*c_ratio) )

        # Adding absolute values of wave speeds to wave_speeds
        np.append(wave_speeds, np.abs(S_L))

        # Setting W_0 (if possible)
        if ( u_star > 0 ):
            if ( S_L > 0 ):
                W_0 = W_L
            else:
                W_0[0], W_0[1], W_0[2] = tuple(rho_star_L, u_star, p_star)
 

    # Obtaining Exact Solution for Right of Contact Discontinuity  
    if ( p_star <= p_R ):
        # Right wave is an expansion wave 
        rho_star_R = rho_R * pow(p_star/p_R, 1/c_ratio)
        S_HR = u_R + a_R
        S_TR = u_star + m.sqrt(c_ratio * p_star / rho_star_R)

        # Adding absolute values of wave speeds to wave_speeds
        np.append(wave_speeds, np.abs([S_HR, S_TR])
        
        # Setting W_0 (if possible)
        if ( u_star < 0 ):
            if ( S_HR < 0 ):
                W_0 = W_R
            elif ( S_TR > 0 ):
                W_0[0], W_0[1], W_0[2] = tuple(rho_star_R, u_star, p_star)
            else:
                W_0[0] = p_R * pow(( (2/(c_ratio+1)) 
                                - (c_ratio-1)/(a_R*(c_ratio+1)) * (u_R) ), 2*c_ratio/(c_ratio-1) )
                W_0[1] = 2/(c_ratio+1) * ( -a_R + u_R * (c_ratio-1)/2 )
                W_0[2] = rho_R * pow(( (2/(c_ratio+1)) 
                                - (c_ratio-1)/(a_R*(c_ratio+1)) * (u_R) ), 2/(c_ratio-1) )
    else:
        # Right wave is a shock wave
        rho_star_R = rho_R * ( (p_star/p_R) 
                        + (c_ratio-1)/(c_ratio+1) ) / ( (c_ratio-1)/(c_ratio+1)*(p_star/p_R) + 1 )
        S_R = u_R + a_R * m.sqrt( (c_ratio+1)/(2*c_ratio) * p_star/p_R + (c_ratio-1)/(2*c_ratio) )

        # Adding absolute values of wave speeds to wave_speeds
        np.append(wave_speeds, np.abs(S_R)
        
        # Setting W_0 (if possible)
        if ( u_star < 0 ):
            if ( S_R < 0 ):
                W_0 = W_R
            else:
                W_0[0], W_0[1], W_0[2] = tuple(rho_star_R, u_star, p_star)


    return W_0, np.max(wave_speeds)



def godunov_method_sys_nonlinear_hyperbolic():
    """
        Function that uses Godunov's Method's second version to obtain a solution of the 
        system of nonlinear, hyperbolic equations over the temporal space
        
        User defines the length of the spatial domain of interest and the number of data
        points to sample data for 
        
        Plots are produced at the end 
    """

    # Specific Heat Ratio 
    c_ratio = 1.4

    # Obtaining user's input
    length = float(input("\nEnter the length of the data line: ")) 
    num_data_pts = int(input("Enter the number of data (the precision to which) the 
                                solution should be obtained for: "))
    dx = length / num_data_pts
    time_period = float(input("\nEnter the time period for which the solution should
                                be obtained"))
    cfl = float(input("CFL coefficient: "))
    print("\nUser should note that the value incremental value for each time step will 
                                be obtained by applying the CFL stability condition."))

    # Creating Primitive Variable Vector, W
    W = np.zeros((3, num_data_pts + 2))
    
    # Initialising the values of W 
    W = initialisation(W, dx)
 
    # Finding how often to plot
    plot_it = int(input("\nNumber of time increments using Godunov's Method should plots
                            be provided: "))


    # GODUNOV'S METHOD
    flux = np.zeros(num_data_pts + 1)
    max_wave_speeds = np.zeros(0)
    X = np.arange(0, length, dx)
    
    t = 0
    count = 1
    while ( t <= time_period ):
        for flx, i in enumerate(flux):
            W_0, max_wave_speed = reimann_problem(W[3,i], W[3,i+1])

            # F(W) = A(W) * W
            flx = np.dot(np.array([[W_0[1], W_0[0], 0],
                                    [0, W_0[1], 1/W_0[0]],
                                    [0, c_ratio*W_0[2], W_0[1]]), W_0)

            # Appending the max wave speed in each Reimann Problem (RP) to max_wave_speeds
            np.append(max_wave_speeds, max_wave_speed)
        
        # Godunov Update Step
        W[3,1:-1] += (1/np.max(max_wave_speeds)) * ( flux[:-1] - flux[1:] ) 
        W[:,0], W[:,-1] = W[:,1], W[:,-2]


        # Incrementing the time 
        t += dx / np.max(max_wave_speeds)

        # ADD plotting feature 
        if ( count % plot_it == 0 ):
            plt.figure()
            
            plt.subplot(311)
            plt.plot(X, W[0,1:-1], 'bo')
            plt.title("Density at t = ",t)
            plt.ylabel("Density")
            plt.xlabel("Position")

            plt.subplot(312)
            plt.plot(X, W[1,1:-1], 'bo')
            plt.title("Velocity at t = ",t)
            plt.ylabel("Velocity")
            plt.xlabel("Position")
            
            plt.subplot(313)
            plt.plot(X, W[2,1:-1], 'bo')
            plt.title("Pressure at t = ",t)
            plt.ylabel("Pressure")
            plt.xlabel("Position")

            plt.tight_layout()
            plt.show()
            
        
        # Incrementing while-loop iteration count 
        count += 1

    

godunov_method_sys_nonlinear_hyperbolic()








 
