import numpy as np 
import math as m
import matplotlib.pyplot as plt
import compressed_erin as cmprin
import slope_limiters as slp
import hllc_reimann_solver as hllc






def RP(U_ebev, M, **kwargs):
    """
        Function that obtains the solution of the RP using the choice Reimann Solver. Calculates
        the intercell flux
        For now, this function just directly uses the HLLC reimann solver. 

        Args:
            U_ebev: ndarray(M+2,2,3); the boundary extrapolated values of the conserved elements
            M: int; number of mesh points
            ...
            **rs: string; reimann solver to be used
            ...

        Returns:
            intercell_flux: ndarray(M+1,3); the flux of the solution to RP at the origin
    """
    # Splitting the evolved BEVs into separate vectors for right & left EBEV for each cell
    ebev_L = U_ebev[:,0]
    ebev_R = U_ebev[:,1]

    # Forming the flux vector
    intercell_flux = np.empty(ebev_R.shape)
    for i in range(M+1):
        intercell_flux[i] = hllc.hllc_1D(ebev_L[i], ebev_R[i])

    return intercell_flux


def stepping_time(U, dt, dx, M):
    """ 
        Function that calculates the values of the conserved variables at the next time step

        Args:
            U: ndarray(M+2,3); the conserved elements 2D array for the time step 
            dt: float; size of the time step 
            dx: float; mesh size
            M: int; number of mesh points

        Returns:
            time_evolved_U: ndarray(M+2,3); the conserved variables for the next time step
    """
    # Obtain the evolved boundary extrapolated values
    U_ebev = bev(U, dt, dx, M)

    # Solving the RPs for the intercell flux
    intercell_flux = RP(U_ebev, M)

    # Calculating conserved variable vector for the next time step
    time_evolved_U = np.empty(U.shape)
    time_evolved_U = U + dt/dx * (intercell_flux[:-1] - intercell_flux[1:])

    # Applying boundary conditions
    time_evolved_U = boundary_conditions(time_evolved_U)

    return time_evolved_U



def main():
    U_bev = np.array([[[0.2,0.1,0.1],
              [0.1,0.1,0.1]],
             [[0.1,0.2,0.2],
              [0.2,0.2,0.2]]])

    U = np.array([[0.1,0.1,0.1],
                  [0.25,0.25,0.25],
                  [0.33,0.33,0.33],
                  [0.42,0.42,0.42],
                  [0.57,0.57,0.57]])

    U_ebev = np.zeros(U_bev.shape)

#    F = np.empty(U_bev.shape)
#    for i, val_1 in enumerate(U_bev):
#        for j, val_2 in enumerate(val_1):
#            F[i,j] = cmprin.flux_euler_1D(val_2)
#    print(F.shape)
    

    # Calculating evolved boundary extrapolated values
#    U_ebev[:,0,:] = U_bev[:,0,:] + 0.5 * (F[:,0,:] - F[:,1,:]) 
#    U_ebev[:,1,:] = U_bev[:,1,:] + 0.5 * (F[:,0,:] - F[:,1,:])

#    print(U_ebev)


    ebev_L = U_bev[:,0]
    ebev_R = U_bev[:,1] 
    print(type(ebev_R))
    print(ebev_R.shape)



if __name__ == "__main__":
    main()


