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



def basic():
    # TODO
    raise NotImplementedError
    return



def xi_R(r, beta, omega):
    """
        Function that calculates the value of xi_R that is required in all the slope
        limters. 

        Args:
            r: float; upwind ratio

        Returns:
            val: float; the calculated value of xi_R
    """

    val = 2 * beta / (1 - omega + (1+omega)*r)
    return val



def xi_L(r):
    # TODO
    raise NotImplementedError
    return



def SUPERBEE(r, beta=1, omega=0):
    """SUPERBEE Slope Limiter 

        Args:
            r: ndarray(3,); upwind ratio
            beta: float; value of beta required for the calculation 
            omega: float; value of omega that is used in the calculation of the simple slope

        Returns:
            xi_SB: float; the SUPERBEE slope limiter
    """

    xi_SB = np.zeros(3)

    for i, val in enumerate(r):
        if ( val <= 0 ):
            xi_SB[i] = 0
        elif ( val <= 0.5 ):
            xi_SB[i] = 2*val
        elif ( val <= 1 ):
            xi_SB[i] = 1
        else:
            xi_SB[i] = min(val, xi_R(val, beta, omega), 2)

    return xi_SB



def van_Leer():
    # TODO
    raise NotImplementedError
    return



def van_Albada():
    # TODO
    raise NotImplementedError
    return



def MINBEE(r, beta=1, omega=0):
    """MINBEE Slope Limiter 

        Args:
            r: ndarray(3,); upwind ratio
            beta: float; value of beta required for the calculation 
            omega: float; value of omega that is used in the calculation of the simple slope

        Returns:
            xi_MB: ndarray(3,); the SUPERBEE slope limiter
    """
   
    xi_MB = np.zeros(3)

    for i, val in enumerate(r):
        if ( val <= 0 ):
            xi_MB[i] = 0
        elif ( val <= 1 ):
            xi_MB[i] = val
        else:
            xi_MB[i] = min(1, xi_R(val, beta, omega))

    return xi_MB



def main():
    # TODO
    raise NotImplementedError
    return



if __name__ == "__main__":
    main()

