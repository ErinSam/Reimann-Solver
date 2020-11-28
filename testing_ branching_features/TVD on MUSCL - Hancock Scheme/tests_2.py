"""
    This testing script is a space to test the reimann solvers and make the supplementary scripts
    required for TVDonMHM work well.
"""


import numpy as np 
import math as m
import compressed_erin as cmprin
from simple_colors import *
import hllc_reimann_solver as hllc








def main():
    # Initialising prim var vectors for testing 
    W_L = np.array([1.0,
                    0.0,
                    1.0])
    W_R = np.array([0.125,
                    0.0,
                    0.1]) 

    # Converting them into the conserved form 
    U_L, U_R = cmprin.prim_consv_1D(W_L), cmprin.prim_consv_1D(W_R)

    # TESTING HLLC 
    F = hllc.hllc_1D(U_L, U_R)
    print(F)

    # TESTING ANRS
#    W_0 = anrs(W_L, W_R)
#    print(W_0)
    print(green("Made it to the end"))








if __name__ == "__main__":
    main()

