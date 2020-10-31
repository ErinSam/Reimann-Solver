"""
    Python Implmentation of an Exact Reimann Solver 

    Version 3: 
        Removing the mistake in the iteration procedure. Previous version started iterating for conditions
        left (and right) of the contact discontinuity at x/t = 0, rather than at x/t = u*, which is the 
        wave speed for the contact discontinuity.
"""
# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electrical & Electronics Minor 



import numpy as np 
import math as m 
import matplotlib.pyplot as plt
from iterative_scheme_for_pressure import *


def exact_reimann_solver():
    
    # INITIALISING VARIABLES
    # Specific Heat Ratio 
    c_ratio = 1.4

    # Denoting the driver region with L and the driven region with R 
    # Obtaining initial conditons for pressure, density and velocity of the driver 
    # and driven region
    p_L = float(input("Enter the value of pressure in the driver region (L): "))
    p_R = float(input("Enter the value of pressure in the driven region (R): ")) 
    rho_L = float(input("Enter the value of density in the driver region (L): "))
    rho_R = float(input("Enter the value of density in the driven region (R): ")) 
    u_L = float(input("Enter the value of velocity in the driver region (L): ")) 
    u_R = float(input("Enter the value of velocity in the driven region (R): "))

    # Obtaining the number of data points and the time step user wants the exact solution for
    # num_data_points = int(input("\nEnter the number of data points you want values to be obtained for: "))
    dx = float(input("\nEnter the value of dx (consider 'dx' as a measure of the accuracy you want the solution to have): "))
    t = float(input("Enter the time step that you want the exact solution for: "))
    
    # Speed of sound in the driver and driven regions
    a_L = m.sqrt(c_ratio * p_L/rho_L)
    a_R = m.sqrt(c_ratio * p_R/rho_R)

    # Initialising the values of the data-dependent constants A_L, B_L, A_R, B_R
    A_L = 2 / ((c_ratio+1)*rho_L)
    B_L = p_L * (c_ratio-1)/(c_ratio+1)
    A_R = 2 / ((c_ratio+1)*rho_R)
    B_R = p_R * (c_ratio-1)/(c_ratio+1)


    # Obtaining the Pressure and Velocity of the Star Region 
    p_star, u_star = star_region_properties(p_L=p_L, p_R=p_R, rho_L=rho_L, rho_R=rho_R, u_L=u_L, u_R=u_R, a_L=a_L, a_R=a_R, A_L=A_L, A_R=A_R, B_L=B_L, B_R=B_R, c_ratio=c_ratio)


    # Primitive Variable Matrix    
    # W[0,] represents pressure
    # W[1,] represents velocity 
    # W[2,] represents density 
    
    # Using two matrices W_L and W_R
    # where W_L is used to obtain and store the values of the primitive variables to 
    # the left of contact discontinuity. Iterations will move from the contact discontinuity 
    # to 50 data points into the left region comprising of initial driver conditions 
    # where W_R is used to obtain and store the value of primitive variables to the right of 
    # the contact discontinuity. Iterations will move from the contact discontinuity 
    # to 50 data points into the right region comprising of initial driven condtions  
    W_L = np.zeros((3, 10000)) 
    W_R = np.zeros((3, 10000))
    
 
    # Obtaining Exact Solution for Left of Contact Discontinuity
    count = 0
    L_count = 0
    self_sim_param = u_star
    if ( p_star <= p_L ):
        # Left wave is an expansion wave 
        print("\nWave to the Left of the Contact Discontinuity is an Expansion Wave")
        rho_star_L = rho_L * pow(p_star/p_L, 1/c_ratio)
        S_HL = u_L - a_L
        S_TL = u_star - m.sqrt(c_ratio * p_star / rho_star_L)
        print("Left Shock Head Speed: ", S_HL)
        print("Left Shock Tail Speed: ", S_TL)

        for count in range(10000):
            if ( self_sim_param <= S_HL ):
                W_L[0, count] = p_L
                W_L[1, count] = u_L
                W_L[2, count] = rho_L
                L_count += 1
                if ( L_count > 50 ): break
            elif ( self_sim_param <= S_TL ):
                W_L[0, count] = p_L * pow(( (2/(c_ratio+1)) + (c_ratio-1)/(a_L*(c_ratio+1)) * (u_L - self_sim_param) ), 2*c_ratio/(c_ratio-1) )
                W_L[1, count] = 2/(c_ratio+1) * ( a_L + u_L * (c_ratio-1)/2 + self_sim_param )
                W_L[2, count] = rho_L * pow(( (2/(c_ratio+1)) + (c_ratio-1)/(a_L*(c_ratio+1)) * (u_L - self_sim_param) ), 2/(c_ratio-1) )
            else: 
                W_L[0, count] = p_star
                W_L[1, count] = u_star
                W_L[2, count] = rho_star_L
        
            count += 1
            self_sim_param -= dx/t

    else:
        # Left wave is a shock wave
        print("\nWave to the Left of the Contact Discontinuity is a Shock Wave")
        rho_star_L = rho_L * ( (p_star/p_L) + (c_ratio-1)/(c_ratio+1) ) / ( (c_ratio-1)/(c_ratio+1)*(p_star/p_L) + 1 )
        S_L = u_L - a_L * m.sqrt( (c_ratio+1)/(2*c_ratio) * p_star/p_L + (c_ratio-1)/(2*c_ratio) )
        print("Left Shock Speed: ", S_L)

        for count in range(10000):
            if ( self_sim_param <= S_L):
                W_L[0, count] = p_L
                W_L[1, count] = u_L
                W_L[2, count] = rho_L
                L_count += 1
                if ( L_count > 50 ): break
            else:
                print("preparing right star region in the case of right shock")
                W_L[0, count] = p_star
                W_L[1, count] = u_star
                W_L[2, count] = rho_star_L
        
            count += 1
            self_sim_param -= dx/t

    # Slicing and flipping the matrix. Preparing it for concatenation
    W_L = W_L[:,:count-1]
    W_L = np.flip(W_L, 1)


    # Obtaining Exact Solution for Right of Contact Discontinuity  
    R_count = 0
    count = 0
    self_sim_param = u_star
    if ( p_star <= p_R ):
        # Right wave is an expansion wave 
        print("\nWave to the Right of the Contact Discontinuity is an Expansion Wave")
        rho_star_R = rho_R * pow(p_star/p_R, 1/c_ratio)
        S_HR = u_R + a_R
        S_TR = u_star + m.sqrt(c_ratio * p_star / rho_star_R)
        print("Right Shock Head Speed: ", S_HR)
        print("Right Shock Tail Speed: ", S_TR)

        for count in range(10000):
            if ( self_sim_param <= S_TR ):
                W_R[0, count] = p_star
                W_R[1, count] = u_star
                W_R[2, count] = rho_star_R
            elif ( self_sim_param <= S_HR ):
                W_R[0, count] = p_R * pow(( (2/(c_ratio+1)) - (c_ratio-1)/(a_R*(c_ratio+1)) * (u_R - self_sim_param) ), 2*c_ratio/(c_ratio-1) )
                W_R[1, count] = 2/(c_ratio+1) * ( -a_R + u_R * (c_ratio-1)/2 + self_sim_param )
                W_R[2, count] = rho_R * pow(( (2/(c_ratio+1)) - (c_ratio-1)/(a_R*(c_ratio+1)) * (u_R - self_sim_param) ), 2/(c_ratio-1) )
            else: 
                W_R[0, count] = p_R
                W_R[1, count] = u_R
                W_R[2, count] = rho_R
                R_count += 1
                if ( R_count > 50 ): break
        
            count += 1
            self_sim_param += dx/t

    else:
        # Right wave is a shock wave
        print("\nWave to the Right of the Contact Discontinuity is a Shock Wave")
        rho_star_R = rho_R * ( (p_star/p_R) + (c_ratio-1)/(c_ratio+1) ) / ( (c_ratio-1)/(c_ratio+1)*(p_star/p_R) + 1 )
        S_R = u_R + a_R * m.sqrt( (c_ratio+1)/(2*c_ratio) * p_star/p_R + (c_ratio-1)/(2*c_ratio) )
        print("Right Shock Speed: ", S_R)

        for count in range(10000):
            if ( self_sim_param <= S_R):
                W_R[0, count] = p_star
                W_R[1, count] = u_star
                W_R[2, count] = rho_star_R
            else:
                W_R[0, count] = p_R
                W_R[1, count] = u_R
                W_R[2, count] = rho_R
                R_count += 1
                if ( R_count > 50 ): break
              
            count += 1
            self_sim_param += dx/t

    # Preparing the matrix for concatenation
    W_R = W_R[:,:count-1]

    
    # Concatenating the matrices
    W = np.concatenate((W_L, W_R), axis=1)
    # Making the x-axis array 
    X = np.array(range(np.shape(W)[1]))

    # Displaying
    print("\n\nPressure in the Star Region: ", p_star)
    print("Velocity in the Star Region: ", u_star)
    print("Density in the L-Star Region: ", rho_star_L)
    print("Density in the R-Star Region: ", rho_star_R)
    print("Values of the primitive variables at the data points:")
    print(W)


    # Plotting
    plt.figure()

    plt.subplot(311)
    plt.plot(X, W[0,])
    plt.title("Exact solution for Pressure")
    plt.ylabel("Pressure")
    plt.xlabel("Position")
    
    plt.subplot(312)
    plt.plot(X, W[1,])
    plt.title("Exact solution for Velocity")
    plt.ylabel("Velocity")
    plt.xlabel("Position")

    plt.subplot(313)
    plt.plot(X, W[2,])
    plt.title("Exact solution for Density")
    plt.ylabel("Density")
    plt.xlabel("Position")
   
    plt.tight_layout()
    plt.show()

#######################################################################################################
#######################################################################################################

exact_reimann_solver()


