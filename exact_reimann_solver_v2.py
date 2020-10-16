# Python Implmentation of an Exact Reimann Solver 
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
    dx = float(input("\nEnter the value of dx: "))
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
    W_L = np.zeros((3, 10000)) 
    W_R = np.zeros((3, 10000))
    
 
    # Obtaining Exact Solution for Left of Contact Discontinuity
    count = 0
    L_count = 0
    self_sim_param = 0
    if ( p_star <= p_L ):
        # Left wave is an expansion wave 
        print("Left wave is an expansion wave")
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
                W_L[2, count] = rho_L * pow(( (2/(c_ratio+1)) + (c_ratio-1)/(a_L*(c_ratio+1)) * (u_L - self_sim_param) ), 2*c_ratio/(c_ratio-1) )
            else: 
                W_L[0, count] = p_star
                W_L[1, count] = u_star
                W_L[2, count] = rho_star_L
        
            count += 1
            self_sim_param -= dx/t

    else:
        # Left wave is a shock wave
        print("Left wave is a shock wave")
        rho_star_L = rho_L * ( (p_star/p_L) + (c_ratio-1)/(c_ratio+1) ) / ( (c_ratio-1)/(c_ratio+1)*(p_star/p_L) + 1 )
        S_L = u_L + a_L * m.sqrt( (c_ratio+1)/(2*c_ratio) * p_star/p_L + (c_ratio-1)/(2*c_ratio) )
        print("Left Shock Speed: ", S_L)

        for count in range(10000):
            if ( self_sim_param <= S_L):
                W_L[0, count] = p_L
                W_L[1, count] = u_L
                W_L[2, count] = rho_L
                L_count += 1
                if ( L_count > 50 ): break
            else:
                W_L[0, count] = p_star
                W_L[1, count] = u_star
                W_L[2, count] = rho_star_L
        
            count += 1
            self_sim_param -= dx/t

    # Slicing and flipping the matrix. Preparing it for concatenation
    print("Number of iteration for L: ", count, L_count)
    print("\n\nW_L looks before splicing like:")
    print(W_L)
    W_L = W_L[:,:count-1]
    W_L = np.flip(W_L, 1)
    print("\n\nW_L looks like:")
    print(W_L)


    # Obtaining Exact Solution for Right of Contact Discontinuity  
    R_count = 0
    count = 0
    self_sim_param = 0
    if ( p_star <= p_R ):
        # Right wave is an expansion wave 
        print("Right wave is an expansion wave")
        rho_star_R = rho_R * pow(p_star/p_R, 1/c_ratio)
        S_HR = u_R + a_R
        S_TR = u_star + m.sqrt(c_ratio * p_star / rho_star_R)
        print("Right Shock Head Speed: ", S_HR)
        print("Right Shock Tail Speed: ", S_TR)

        for count in range(10000):
            if ( self_sim_param <= S_HR ):
                W_R[0, count] = p_star
                W_R[1, count] = u_star
                W_R[2, count] = rho_star_R
            elif ( self_sim_param <= S_TR ):
                W_R[0, count] = p_R * pow(( (2/(c_ratio+1)) + (c_ratio-1)/(a_R*(c_ratio+1)) * (u_R - self_sim_param) ), 2*c_ratio/(c_ratio-1) )
                W_R[1, count] = 2/(c_ratio+1) * ( -a_R + u_R * (c_ratio-1)/2 + self_sim_param )
                W_R[2, count] = rho_R * pow(( (2/(c_ratio+1)) + (c_ratio-1)/(a_R*(c_ratio+1)) * (u_R - self_sim_param) ), 2*c_ratio/(c_ratio-1) )
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
        print("Right wave is a shock wave")
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
    print("Number of iteration for R: ", count, R_count)
    print("\n\nW_R looks before splicing like:")
    print(W_R)
    W_R = W_R[:,:count-1]
    print("\n\nW_R looks like:")
    print(W_R)

    
    # Concatenating the matrices
    W = np.concatenate((W_L, W_R), axis=1)
    # Making the x-axis array ... is this necessary?
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
    plt.ylabel("Pressure")
    
    plt.subplot(312)
    plt.plot(X, W[1,])
    plt.ylabel("Velocity")

    plt.subplot(313)
    plt.plot(X, W[2,])
    plt.ylabel("Density")
   
    plt.show()

#######################################################################################################
#######################################################################################################

exact_reimann_solver()


