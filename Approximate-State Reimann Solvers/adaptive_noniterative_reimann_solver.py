"""
    Python script for the Adaptive Noniterative Reimann Solver (ANRS), an approximate-
    state reimann solver. 

    Hybrid schemes work on the basis of the use of a simple Reimann Solver in regions 
    of smooth flow (almost similar conditions) and near isolated contacts and shear 
    waves, and a more sophisticated Reimann Solver elsewhere, in an adaptive fashion, 
    which is defined by the user's choice of a switchin parameter
"""

# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electronics & Electrical Engg. Minor


import numpy as np 
import math as m 


c_ratio = 1.4

def pvrs(W_L, W_R, **kwargs):
    """ PRIMITIVE VARIABLE REIMANN SOLVER 

        Primitive variable reimann solver that allows the approximation to be made as a
        default using averages or using characteristic equations. 

        Args:
            W_L: ndarray(3,); Primitve variable vector for the left region (driver region)
            W_R: ndarray(3,); Primitive variable vector for the region region (driven region)
            kwargs: 
                'charac_eqns': Bool; Obtains approximate solutions for the star values using 
                    characterisitic equations
                default: Obtains approximate solutions for the star values using averages

        Returns:
            p_star: float; Star region pressure 
            u_star: float; Star region velocity 
            rho_L_star: float; Left star region density
            rho_R_star: float; Right star region density
    """
    
    rho_L, u_L, p_L = W_L[0], W_L[1], W_L[2]
    rho_R, u_R, p_R = W_R[0], W_R[1], W_R[2]
    a_L, a_R = m.sqrt(c_ratio * p_L / rho_L), m.sqrt(c_ratio * p_R / rho_R)

    if ( kwargs['charac_eqns'] == True ):
        C_L = rho_L * a_L
        C_R = rho_R * a_R

        # Calculating Star Region Properties
        p_star = 1/(C_L + C_R) * ( C_R*p_L + C_L*p_R + C_L*C_R*(u_L-u_R) )
        u_star = 1/(C_L + C_R) * ( C_L*u_L + C_R*u_R + (p_L - p_R) )
        rho_L_star = rho_L + (p_star - p_L) / pow(a_L,2)
        rho_R_star = rho_R + (p_star - p_R) / pow(a_R,2) 
    
    else:
        rho_avg = 0.5 * (rho_L + rho_R)
        a_avg = 0.5 * (a_L + a_R)

        # Calculating Star Region Properties
        p_star = 0.5 * (p_L + p_R) + 0.5 * (u_L - u_R) * (rho_avg * a_avg)
        u_star = 0.5 * (u_L + u_R) + 0.5 * (p_L - p_R) * (rho_avg * a_avg)
        rho_L_star = rho_L + ( u_L - u_star ) * (rho_avg/a_avg) 
        rho_R_star = rho_R + ( u_star - u_R ) * (rho_avg/a_avg)

    return p_star, u_star, rho_L_star, rho_R_star 
        


def trrs(W_L, W_R):
    """ TWO RAREFACTION REIMANN SOLVER 

        Obtaining the star region values in the case that both the non-linear waves are
        rarefaction waves. 
        Star region density & velocity are obtained using exact relations. 

        Args:
            W_L: ndarray(3,); Primitve variable vector for the left region (driver region)
            W_R: ndarray(3,); Primitive variable vector for the region region (driven region)

        Returns:
            p_star: float; Star region pressure 
            u_star: float; Star region velocity 
            rho_L_star: float; Left star region density
            rho_R_star: float; Right star region density
    """

    rho_L, u_L, p_L = W_L[0], W_L[1], W_L[2]
    rho_R, u_R, p_R = W_R[0], W_R[1], W_R[2]
    a_L, a_R = m.sqrt(c_ratio * p_L / rho_L), m.sqrt(c_ratio * p_R / rho_R)

    z = (c_ratio - 1) / (2*c_ratio)

    def f_k(P, W):
        """ This is a generalized f_L & f_R functions in the pressure function
        """
        rho, u, p = W[0], W[1], W[2]
        a = m.sqrt(c_ratio * p / rho)
        A = 2 / ( (c_ratio+1)*rho )
        B = p * (c_ratio - 1)/(c_ratio+1)
        
        if ( P > p ):
            return ( (P-p) * pow(A/(P+B), 0.5) )
        else:
            return ( 2*a/(c_ratio-1) * (pow(P/p,z)-1) )

    # Calculating Star Region Properties
    p_star = m.sqrt( (a_L + a_R - (c_ratio-1)/2 * (u_R-u_L)) \
                / ( a_L/pow(p_L,z) + a_R/pow(p_R,z) ), 1/z)
    u_star = 0.5 * (u_L + u_R) + 0.5 ( f_k(p_star, W_R) - f_k(p_star, W_L) )
    rho_L_star = rho_L * pow(p_star/p_L, 1/c_ratio)
    rho_R_star = rho_R + pow(p_star/p_R, 1/c_ratio) 

    return p_star, u_star, rho_L_star, rho_R_star 



def tsrs(W_L, W_R, p_0):
    """ TWO SHOCK REIMANN SOLVER 

        Obtaining the star region values in the case that both the non-linear waves are
        shock waves. 

        Args:
            W_L: ndarray(3,); Primitve variable vector for the left region (driver region)
            W_R: ndarray(3,); Primitive variable vector for the region region (driven region)
            p_0: float; Estimate for the solution of pressure: max(0, p_pvrs)

        Returns:
            p_star: float; Star region pressure 
            u_star: float; Star region velocity 
            rho_L_star: float; Left star region density
            rho_R_star: float; Right star region density
    """

    rho_L, u_L, p_L = W_L[0], W_L[1], W_L[2]
    rho_R, u_R, p_R = W_R[0], W_R[1], W_R[2]
    a_L, a_R = m.sqrt(c_ratio * p_L / rho_L), m.sqrt(c_ratio * p_R / rho_R)

    z = (c_ratio - 1) / (2*c_ratio)
    weird_ratio = (c_ratio - 1)/(c_ratio + 1)

    def g_k(P, W):
        """ This is another generalized function that is present as a term in the pressure 
            function
        """
        rho, u, p = W[0], W[1], W[2]
        a = m.sqrt(c_ratio * p / rho)
        A = 2 / ( (c_ratio+1)*rho )
        B = p * (c_ratio - 1)/(c_ratio+1)

        return pow(A/(P + B), 0.5)

    g_L = g_k(p_0, W_L)
    g_R = g_k(p_0, W_R)

    # Calculating Star Region Properties
    p_star = ( g_L*p_L + g_R*p_R - (u_R - u_L) ) / ( g_L + g_R )
    u_star = 0.5 * (u_L + u_R) + 0.5*( (p_star-p_R)*g_R - (p_star-p_L)*g_L )
    rho_L_star = ( p_star/p_L + weird_ratio ) / ( weird_ratio * p_star/p_L + 1)
    rho_R_star = ( p_star/p_R + weird_ratio ) / ( weird_ratio * p_star/p_R + 1)

    return p_star, u_star, rho_L_star, rho_R_star 




def anrs(W_L, W_R, **kwargs):
    """ ADAPTIVE NON-ITERATIVE REIMANN SOLVER

        Function that solves a Reimann Problem using an Adaptive Noniterative Reimann
        Solver. A combination of the PVRS scheme as the cheap component and the TRRS 
        and TSRS solvers to provide the robust component of the adaptive scheme. 

        The PVRS scheme is used in the case that:
            (i) p_max / p_min < Switching Parameter
            (ii) p_min < Star Region Pressure < p_max

        Args:
            W_L: ndarray(3,); Primitve variable vector for the left region (driver region)
            W_R: ndarray(3,); Primitive variable vector for the region region (driven region)
            kwargs:
                'pvrs':

        Returns:
    """

    
    # Assigning value for the Switching Parameter
    switching_param = 2 

    
