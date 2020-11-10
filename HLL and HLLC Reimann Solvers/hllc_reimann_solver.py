"""
    Python script that implements the HLLC Reimann Solver, a variation of the HLL Reimann Solver
    that was developed by Toro and contributors. 

    For now, this script only implmenets the HLLC RS to solve RPs for 1D Euler equations 
"""

# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electrical & Electronics Engg. Minor



import numpy as np 
import math as m
from adaptive_noniterative_reimann_solver import pvrs



def hllc(W_L, W_R):
    """ HLLC REIMANN SOLVER 

        This function solves a Reimann problem given left and right region properties 
        and specifically obtains the solution at the origin, that is, the point where
        the similarity parameter is zero. 

        HLLC RS can be used within Godunov Schemes to find the intercell numericall flux

        Args:
            W_L: ndarray(3,); Primitve variable vector for the left region (driver region)
            W_R: ndarray(3,); Primitive variable vector for the region region (driven region)

        Returns:
            F_0: ndarray(3,); Intercell Numerical Flux Vector by HLLC @ x/t = 0 
    """
    
    # Heat Constant Ratio
    c_ratio = 1.4

    # Obtaining parameters from primitive variable vectors 
    rho_L, u_L, p_L = W_L[0], W_L[1], W_L[2]
    rho_R, u_R, p_R = W_R[0], W_R[1], W_R[2]
    a_L, a_R = m.sqrt(c_ratio * p_L / rho_L), m.sqrt(c_ratio * p_R / rho_R)
   
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
    # TODO     
    def U_k():
        # TODO

    def F_k_star():
        # TODO

    # Location of the required intercell flux 
    # TODO


    # Return Statement 
    # TODO
