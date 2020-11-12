"""
    Python script that implements the Reimann Solver of Roe.

    So far, this script is only capable of implementing the Roe Reimann Solver to solve 
    RPs for x-split 3D Euler Equations.  
"""

# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electronics & Electrical Engg. Minor



import numpy as np 
import math as m



def roe_pike_reimann_solver(U_L, U_R):
    """ ROE'S REIMANN SOLVER 

        This function uses the reimann solver of roe to solve Reimann Problems for x-split 
        3D Euler Equations. 
        This version of Roe's Reimann solver was introduced by Roe & Pike. This is does not 
        involve the entropy fix. 

        Args:
            U_L: ndarray(5,); vector of L conserved variable for the x-split 3D Euler equations
            U_R: ndarray(5,); vector of R conserved variable for the x-split 3D Euler equations

        Returns: 
            intercell_flux: ndarray(5,); vector of intercell numerical flux that is calculated 
                                            by the Roe Reimann Solver
    """

    # Heat Capacity Ratio 
    c_ratio = 1.4

    # Obtain the flow field variables from the conserved vectors 
    rho_L, u_L, v_L, w_L, E_L = U_L[0], U_L[1]/U_L[0], U_L[2]/U_L[0], U_L[3]/U_L[0], U_L[4] 
    p_L = (c_ratio-1) * ( E_L - 0.5 * rho_L * (pow(u_L,2) + pow(v_L,2) + pow(w_L,2) ) )
    H_L = (E_L + p_L) / rho_L

    rho_R, u_R, v_R, w_R, E_R = U_R[0], U_R[1]/U_R[0], U_R[2]/U_R[0], U_R[3]/U_R[0], U_R[4] 
    p_R = (c_ratio-1) * ( E_R - 0.5 * rho_R * (pow(u_R,2) + pow(v_R,2) + pow(w_R,2) ) )
    H_R = (E_R + p_R) / rho_R


    # Computing the Roe - Averaged Variables 
    rho_avg = m.sqrt(rho_L * rho_R)
    u_avg = ( m.sqrt(rho_L)*u_L + m.sqrt(rho_R)*u_R ) / ( m.sqrt(rho_L) + m.sqrt(rho_R) )
    v_avg = ( m.sqrt(rho_L)*v_L + m.sqrt(rho_R)*v_R ) / ( m.sqrt(rho_L) + m.sqrt(rho_R) )
    w_avg = ( m.sqrt(rho_L)*w_L + m.sqrt(rho_R)*w_R ) / ( m.sqrt(rho_L) + m.sqrt(rho_R) )
    H_avg = ( m.sqrt(rho_L)*H_L + m.sqrt(rho_R)*H_R ) / ( m.sqrt(rho_L) + m.sqrt(rho_R) )
    a_avg = m.sqrt( (c_ratio-1) * (H_avg - 0.5 * (pow(u_avg,2) + pow(v_avg,2) + pow(w_avg,2)) 


    # Computing the averaged eigenvalues using the analytical expressions 
    lam_1_avg,lam_2_avg,lam_3_avg,lam_4_avg,lam_5_avg = (u_avg-a_avg), u_avg, u_avg, u_avg, (u_avg+a_avg)
    
    # Creating an averaged lambda vector for ease of calculation
    lam_avg_vec = np.array([lam_1_avg,
                            lam_2_avg,
                            lam_3_avg,
                            lam_4_avg,
                            lam_5_avg])

    # Computing the averaged right eigenvectors using the analytical expressions 
    K_1_avg = np.array([1,
                        u_avg - a_avg,
                        v_avg,
                        w_avg,
                        H_avg - u_avg*a_avg])
    K_2_avg = np.array([1,
                        u_avg,
                        v_avg,
                        w_avg,
                        0.5 * (pow(u_avg,2) + pow(v_avg,2) + pow(w_avg,2))])
    K_3_avg = np.array([0,
                        0,
                        1,
                        0,
                        v_avg])
    K_4_avg = np.array([0,
                        0,
                        0,
                        1,
                        w_avg])
    K_5_avg = np.array([1,
                        u_avg + a_avg,
                        v_avg,
                        w_avg,
                        H_avg + u_avg*a_avg])

    # Creating an averaged right eigenvector vector for ease of calculation
    K_avg_vec = np.array([K_1_avg,
                            K_2_avg,
                            K_3_avg,
                            K_4_avg,
                            K_5_avg]) 

    # Computing the averaged wavestrengths using the analytical expressions 
    du, dv, dw, dp, drho = u_R - u_L, v_R - v_L, w_R - w_L, p_R - p_L, rho_R - rho_L

    alpha_1_avg = 1/(2*pow(a_avg,2)) * (dp - rho_avg*a_avg*du)
    alpha_2_avg = drho - dp / pow(a_avg,2)
    alpha_3_avg = rho_avg * dv
    alpha_4_avg = rho_avg * dw
    alpha_5_avg = 1/(2*pow(a_avg,2)) * (dp + rho_avg*a_avg*du) 
    
    # Creating an averaged wavestrength vector for ease of calculation
    alpha_avg_vec = np.array([alpha_1_avg,
                                alpha_2_avg,
                                alpha_3_avg,
                                alpha_4_avg,
                                alpha_5_avg])
    
    # Computing the intercell numerical flux 
    # Left Flux, F(U_L)
    F_L = np.array([rho_L * u_L,
                    rho_L * pow(u_L,2) + p_L,
                    rho_L * u_L * v_L,
                    rho_L * u_L * w_L,
                    u_L * (E_L + p_L)])

    # Right Flux, F(U_R)
    F_R = np.array([rho_R * u_R,
                    rho_R * pow(u_R,2) + p_R,
                    rho_R * u_R * v_R,
                    rho_R * u_R * w_R,
                    u_R * (E_R + p_R)])

    # Shift Flux 
    flux_sum_vec = np.zeros(5)
    for i, row in enumerate(K_avg_vec):
        flux_sum_vec += (alpha_avg_vec[i] * abs(lam_avg_vec) * row)
        
    # Intercell Flux, F(U_0) 
    F_0 = 0.5 * (F_L + F_R) - 0.5 * flux_sum_vec
    
    
    # Return statement 
    return F_0



def entropy_fix_roe_pike_reimann_solver(U_L, U_R, head_speed, tail_speed, rho_star, direction):
    """ ENTROPY FIX 

        This function implements the Harten-Hyman Entropy Fix that is required to solve the 
        RP within a rarefaction
 
        Args:
            U_L: ndarray(5,); vector of L conserved variable for the x-split 3D Euler equations
            U_R: ndarray(5,); vector of R conserved variable for the x-split 3D Euler equations
            head_speed: float; speed of the head as calculated from the chosen form of star
                                region approximation
            tail_speed: float; speed of the head as calculated from the chosen form of star
                                region approximation 
            rho_star: float; 
            direction: string; 

        Returns: 
            intercell_flux: ndarray(5,); vector of intercell numerical flux that is calculated 
                                            by the Roe Reimann Solver
    """

    # Obtaining the flow-field variables from the conserved variabl vectors
    rho_L, u_L, v_L, w_L, E_L = U_L[0], U_L[1]/U_L[0], U_L[2]/U_L[0], U_L[3]/U_L[0], U_L[4] 
    p_L = (c_ratio-1) * ( E_L - 0.5 * rho_L * (pow(u_L,2) + pow(v_L,2) + pow(w_L,2) ) )
    H_L = (E_L + p_L) / rho_L

    rho_R, u_R, v_R, w_R, E_R = U_R[0], U_R[1]/U_R[0], U_R[2]/U_R[0], U_R[3]/U_R[0], U_R[4] 
    p_R = (c_ratio-1) * ( E_R - 0.5 * rho_R * (pow(u_R,2) + pow(v_R,2) + pow(w_R,2) ) )
    H_R = (E_R + p_R) / rho_R


    # Calculation of the Roe - Averaged Variables
    rho_avg = m.sqrt(rho_L * rho_R)
    u_avg = ( m.sqrt(rho_L)*u_L + m.sqrt(rho_R)*u_R ) / ( m.sqrt(rho_L) + m.sqrt(rho_R) )
    v_avg = ( m.sqrt(rho_L)*v_L + m.sqrt(rho_R)*v_R ) / ( m.sqrt(rho_L) + m.sqrt(rho_R) )
    w_avg = ( m.sqrt(rho_L)*w_L + m.sqrt(rho_R)*w_R ) / ( m.sqrt(rho_L) + m.sqrt(rho_R) )
    H_avg = ( m.sqrt(rho_L)*H_L + m.sqrt(rho_R)*H_R ) / ( m.sqrt(rho_L) + m.sqrt(rho_R) )
    a_avg = m.sqrt( (c_ratio-1) * (H_avg - 0.5 * (pow(u_avg,2) + pow(v_avg,2) + pow(w_avg,2)) 


    # Calculating the required averaged right eigenvector and averaged wavestrength
    du, dv, dw, dp, drho = u_R - u_L, v_R - v_L, w_R - w_L, p_R - p_L, rho_R - rho_L
    if ( direction == "left" ):
        K = np.array([1,
                        u_avg - a_avg,
                        v_avg,
                        w_avg,
                        H_avg - u_avg*a_avg])
        alpha = 1/(2*pow(a_avg,2)) * (dp - rho_avg*a_avg*du)
    
    else:
        K = np.array([1,
                        u_avg + a_avg,
                        v_avg,
                        w_avg,
                        H_avg + u_avg*a_avg])
        alpha = 1/(2*pow(a_avg,2)) * (dp + rho_avg*a_avg*du) 


    # Calulation of the Roe-Averaged wave speed between the head and tail of the shock wave 
    if ( direction == "left" ):
        wave_avg = ( m.sqrt(rho_L)*tail_speed + m.sqrt(rho_star)*head_speed ) \
                    / ( m.sqrt(rho_L) + m.sqrt(rho_star) ) 
    else:
        wave_avg = ( m.sqrt(rho_star)*tail_speed + m.sqrt(rho_R)*head_speed ) \
                    / ( m.sqrt(rho_star) + m.sqrt(rho_R) ) 


    # Calculation of the new wave speed
    if ( direction == "left" ):
        new_wave_speed = tail_speed * ( head_speed - wave_avg ) / ( head_speed - tail_speed )
    else:
        new_wave_speed = head_speed * ( wave_avg - tail_speed ) / ( head_speed - tail_speed )


    # Computing the intercell numerical flux 
    # Left Flux, F(U_L)
    F_L = np.array([rho_L * u_L,
                    rho_L * pow(u_L,2) + p_L,
                    rho_L * u_L * v_L,
                    rho_L * u_L * w_L,
                    u_L * (E_L + p_L)])

    # Right Flux, F(U_R)
    F_R = np.array([rho_R * u_R,
                    rho_R * pow(u_R,2) + p_R,
                    rho_R * u_R * v_R,
                    rho_R * u_R * w_R,
                    u_R * (E_R + p_R)])
    

    # Obtaining the intercell numerical flux 
    if ( direction == "left" ):
        intercell_flux = F_L + new_wave_speed * alpha * K
    else:
        intercell_flux = F_R - new_wave_speed * alpha * K 


    # Return statement
    return intercell_flux








def roe_pike_reimann_solver_main(U_L, U_R):
    """ ROE PIKE REIMANN SOLVER 

        This funtion implements the Roe-Pike method to solve a RP given the L and R region 
        properties. It includes the entropy fix that is used to check whether the solution 
        for the intercell numerical flux is present within a rarefaction wave. 

        Args:
            U_L: ndarray(5,); vector of L conserved variable for the x-split 3D Euler equations
            U_R: ndarray(5,); vector of R conserved variable for the x-split 3D Euler equations

        Returns: 
            intercell_flux: ndarray(5,); vector of intercell numerical flux that is calculated 
                                            by the Roe Reimann Solver
    """

    # Heat Capacity Ratio 
    c_ratio = 1.4

    # Obtaining approximation for the star region with ANRS
    # TODO

    # Checking whether entropy fix should be used or not 
    # TODO

    # Calling the right functions 
    # TODO

    # Return statement
    # TODO


 
