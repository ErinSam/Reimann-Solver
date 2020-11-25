"""
    Python script (more of a module(?)) that contains a variety of slope limiters that can be
    used in TVD or high-order schemes

    All slope limiters have been currently implemented to handle the 1D Euler Equations in 
    the conserved form 
        U_t + F(U)_x = 0
"""
# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electronics & Electrical Engg. Minor



import numpy as np 
import math as m



if __name__ == "__main__":
    main()



def basic():
    # TODO



def xi_R(r):
    # TODO




def xi_L(r):
    # TODO


def SUPERBEE(r, xi_R):
    """SUPERBEE Slope Limiter 

        Args:
            r: float; upwind ratio
            xi_R: function; the xi_R function that is required for calculation

        Returns:
            xi_SB: float; the SUPERBEE slope limiter
    """

    if ( r <= 0 ):
        return 0
    elif ( r <= 0.5 ):
        return 2*r
    elif ( r <= 1 ):
        return 1
    else: 
        return min(r, xi_R(r), 2)




def van_Leer():
    # TODO



def van_Albada():
    # TODO



def MINBEE(r, xi_R):
    """MINBEE Slope Limiter 

        Args:
            r: float; upwind ratio
            xi_R: function; the xi_R function that is required for calculation

        Returns:
            xi_SB: float; the SUPERBEE slope limiter
    """
    
    if ( r <= 0 ):
        return 0
    elif ( r <= 1 ):
        return r
    else:
        return min(1, xi_R)



def main():
    # TODO
