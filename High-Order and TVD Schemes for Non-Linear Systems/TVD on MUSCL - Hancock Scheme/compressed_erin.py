"""
    This python script contains functions that kept showing up in scripts that I had to write,
    in one way or another. 

    1. consv_prim_1D(U)
        Function that converts a 1D Euler Equations conserved variables vector into its 
        corresponding primtive variables vector.
        Under the assumption that the gas that we are dealing with obeys the Ideal Gas EOS. 

        Args:
            U: ndarray(3,); conserved variables vector

        Returns:
            W: ndarray(3,); primtive variables vector
"""
# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electronics & Electrical Engg. Minor


import numpy as np 
import math as m



def consv_prim_1D(U):
    """
        Function that converts a 1D Euler Equations conserved variables vector into its 
        corresponding primtive variables vector.
        Under the assumption that the gas that we are dealing with obeys the Ideal Gas EOS. 

        Args:
            U: ndarray(3,); conserved variables vector

        Returns:
            W: ndarray(3,); primtive variables vector
    """
    c_ratio = 1.4

    W = np.array([U[0],
                    U[1]/U[0],
                    (c_ratio-1)*U[0]*(U[2]/U[0] - 0.5*pow(U[1]/U[0],2))])

    return W


