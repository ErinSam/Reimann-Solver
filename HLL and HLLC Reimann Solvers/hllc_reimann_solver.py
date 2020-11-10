"""
    Python script that implements the HLLC Reimann Solver, a variation of the HLL Reimann Solver
    that was developed by Toro and contributors. 

    For now, this script only implmenets the HLLC RS to solve RPs for 1D Euler equations 
"""

# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electrical & Electronics Engg. Minor



import numpy as np 
import math as m
from adaptive_noniterative_reimann_solver import pvrs



def hllc(U_L, U_R):
    """ HLLC REIMANN SOLVER 

        This function solves a Reimann problem given left and right region properties 
        and specifically obtains the solution at the origin, that is, the point where
        the similarity parameter is zero. 

        HLLC RS can be used within Godunov Schemes to find the intercell numericall flux

        Args:
            U_L: ndarray(3,); Flow field variable vector for the left region (driver region)
            U_R: ndarray(3,); Flow field variable vector for the region region (driven region)

        Returns:
            F_0: ndarray(3,); Intercell Numerical Flux Vector by HLLC @ x/t = 0 
    """
    
    # Heat Constant Ratio
    c_ratio = 1.4

    # Obtaining parameters from primitive variable vectors 
    rho_L, u_L, p_L = U_L[0], U_L[1]/U_L[0], ( (c_ratio-1)*U_L[2] - 0.5*pow(U_L[1],2)/U_L[0] )
    rho_R, u_R, p_R = U_R[0], U_R[1]/U_R[0], ( (c_ratio-1)*U_R[2] - 0.5*pow(U_R[1],2)/U_R[0] )
    a_L, a_R = m.sqrt(c_ratio * p_L / rho_L), m.sqrt(c_ratio * p_R / rho_R)
   
    # Initialising primitive variable vectors
    W_L = np.array([rho_L, u_L, p_L]) 
    W_R = np.array([rho_R, u_R, p_R]) 

    # Pressure Estimate 
    p_star, u_star, rho_L_star, rho_R_star = pvrs(W_L, W_R)
    p_star = max(0, p_star)

    # Wave Speed Estimates
    def q_k(p):
        """ Generalized function for k = L, R
            Does a quick calculation on p_L or p_R
        """
        if ( p_star <= p ):
            return 1
        else:
            return pow( (1 + (c_ratio+1)/(2*c_ratio) * (p_star/p - 1)), 0.5)

    S_L = u_L - a_L * q_k(p_L)
    S_R = u_R + a_R * q_k(p_R)
    S_star = ( p_R - p_L + rho_L*u_L*(S_L - u_L) - rho_R*u_R*(S_R - u_R) )\
                / ( rho_L*(S_L - u_L) - rho_R*(S_R - u_R) )


    # HLLC Flux
    def U_k_star(U_k, W_k, S_k):
        """ Generalized macro used to obtain the conserved variable vector for the star region

            Args:
                U_k: ndarray(3,); Conserved variable vector of the k-region 
                W_k: ndarray(3,); Primitive variable vector of the k-region 
                S_k: float; Wave speed of the non-linear wave adjacent to k-region

            Returns:
                U_k_star: ndarray(3,); Conserved variable vector of k-star region 
                
        """
        wave_speed_coeff = (S_k - W_k) / (S_k - S_star)
        star_region_comp_vec_energy_term = U_k[2]/W_k[0] + (S_star - W_k[1]) \
                                            * ( S_star + W_k[2] / (W_k[0]*(S_k - W_k[1])) )
        star_region_comp_vec = np.array([1, S_star, star_region_comp_vec_energy_term])

        star_region_var_vec = W_k[0] * wave_speed_coeff * star_region_comp_vec 
        return star_region_var_vec 

    def F_k_star(U_k, W_k, S_k, F_k):
        """ Generalized macro to calculate k-star region intercell numerical flux

            Args:
                U_k: ndarray(3,); Conserved variable vector of the k-region 
                W_k: ndarray(3,); Primitive variable vector of the k-region 
                S_k: float; Wave speed of the non-linear wave adjacent to k-region
                F_k: ndarray(3,); Intercell numerical flux of the k-region

            Returns:
                star_region_flux: ndarray(3,); k-star region intercell numerical flux 
        """
        D_star = np.array([0, 1, S_star])
        star_region_flux = ( S_star * (S_k * U_k - F_k) \
                            + S_k * (W_k[2] + rho_L*(S_k - W_k[1])*(S_star - W_k[1])) * D_star ) \
                            / (S_k - S_star)
        return star_region_flux

    # Location of the required intercell flux 
    if ( 0 <= S_L ):
        # Flux is calculated from the LEFT region
        f_1 = rho_L * u_L
        f_2 = rho_L * pow(u_L,2) + p_L
        f_3 = u_L * (U_L[2] + p_L)

        F_L = np.array([f_1, f_2, f_3])
        F = F_L
    
    elif ( S_L <= 0 & S_star >= 0 ):
        # Flux is calculated from the LEFT STAR region
        f_1 = rho_L * u_L
        f_2 = rho_L * pow(u_L,2) + p_L
        f_3 = u_L * (U_L[2] + p_L)

        F_L = np.array([f_1, f_2, f_3])
        F = F_L + S_L * (U_k_star(U_L, W_L, S_L) - U_L)

    elif ( S_star <= 0 & S_R >= 0 ):
        # Flux is calculated from the RIGHT STAR region
        f_1 = rho_R * u_R
        f_2 = rho_R * pow(u_R,2) + p_R
        f_3 = u_R * (U_R[2] + p_R)

        F_R = np.array([f_1, f_2, f_3])
        F = F_R + S_R * (U_k_star(U_R, W_R, S_R) - U_R)

    else: 
        # Flux is calculated from the RIGHT region
        f_1 = rho_R * u_R
        f_2 = rho_R * pow(u_R,2) + p_R
        f_3 = u_R * (U_R[2] + p_R)

        F_R = np.array([f_1, f_2, f_3])
        F = F_R


    # Return Statement 
    return F
