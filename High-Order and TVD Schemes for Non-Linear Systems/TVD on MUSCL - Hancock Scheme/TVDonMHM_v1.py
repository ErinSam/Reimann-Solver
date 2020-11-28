"""
    Python script that implements Total Variation Diminishing (TVD) on the MUSCL-Hancock
    scheme to solve 1D Euler Equations
"""
# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electrical & Electronics Engg. Minor



import numpy as np 
import math as m
import matplotlib.pyplot as plt
import slope_limiters as slp
import compressed_erin as cmprin
import hllc_reimann_solver as hllc



def initialisation(U, L, M):
    """
        Function that initialises the initial data line

        Args:
            U: ndarray(M+2,3); the conserved elements 2D array for the time step 
            L: float; length of the data line
            M: int; number of mesh points

        Returns:
            U: ndarray(M+2,3); the conserved elements 2D array for the time step after applying 
                                the initial coniditons to it
    """
    # TODO


    # Defintion of the initial conditions
    raise NotImplementedError
    return



def boundary_conditions(U):
    """
        Function that applies the boundary conditions according to our choice of 
        boundary conditions
        
        Args: 
            U: ndarray(M+2,3); the conserved elements 2D array for the time step 

        Returns: 
            U: ndarray(M+2,3); the conserved elements 2D array for the time step after applying 
                                the boundary coniditons to it
    """
    # Transmissive Boundary Conditions
    U[0], U[-1] = U[1], U[-2]

    # Relfective Boundary Conditions 
    # TODO

    return U



def bev(U, dt, dx, M, omega=0.5):
    """ BOUNDARY EXTRAPOLATED VALUES
        Function that calculates the boundary extrapolated values. Results are used to obtain the
        evolved BEVs and they depend a lot on the choice of the slope limter chosen.

        CAUTION: I have no idea what to take omega as!!

        Args:
            U: ndarray(M+2,3); the conserved elements 2D array for the time step 
            dt: float; size of the time step 
            dx: float; mesh size
            M: int; number of mesh points

        Returns:
    """
    # Creating ndarray for boundary extrapolated values
    U_bev = np.zeros((M+2,2,3))

    # Obtaining the intercell element-wise jump vector (intrcll_delta)
    intercell_data = U[1:] - U[:-1]

    # Calculating upwind ratio
    r = np.empty(((U_bev.shape[0] - 2), 3))
    for i in range(r.shape[0]):
        for j in range(3):
            if ( (intercell_data[i,j] == 0) | (intercell_data[i+1,j] == 0) ):
                r[i,j] = 0
            else:
                r[i,j] = intercell_data[i,j] / intercell_data[i+1,j]

    # Calculating the initial slopes
    delta = 0.5 * (1+omega) * intercell_data[:-1] + 0.5 * (1-omega) * intercell_data[1:]

    # Obtaining slope limiter
    XI = np.zeros(r.shape)
    for i, val in enumerate(r):
        XI[i] = slp.SUPERBEE(val)                       # FIXME

    # Calculcating the limited slopes
    limited_slope = XI * delta
    
    # Obtaining boundary extrapolated values
    U_bev[1:-1,0] = U[1:-1] - limited_slope/2
    U_bev[1:-1,1] = U[1:-1] + limited_slope/2

    # Applying zero slope boundary conditions to the boundary BEVs
    U_bev[0,:] = U[0]
    U_bev[-1,:] = U[-1]

    # Obtaining Evolved BEVs
    U_ebev = evolved_bev(U_bev, dt, dx)

    return U_ebev



def evolved_bev(U_bev, dt, dx):
    """
        Function that caluclates the evolved boundary extrapolated values 
        This results of this function are then used to obtain the intercell flux by solvng the RP

        Args:
            U_bev: ndarray(M+2,2,3); the boundary extrapolated values of the conserved elements
            dt: float; size of the time step 
            dx: float; mesh size

        Returns:
            U_ebev: ndarray(M+2,2,3); evolved boundary extrapolated values of the conserved vars
    """
    # Creatin ndarray for eveolved boundary extrapolated values
    U_ebev = np.zeros(U_bev.shape)

    # Calculating the flux function
    F = np.empty(U_bev.shape)
    for i, val_1 in enumerate(U_bev):
        for j, val_2 in enumerate(val_1):
            F[i,j] = cmprin.flux_euler_1D(val_2)

    # Calculating evolved boundary extrapolated values
    U_ebev[:,0,:] = U_bev[:,0,:] + 0.5*dt/dx * (F[:,0,:] - F[:,1,:])
    U_ebev[:,1,:] = U_bev[:,1,:] + 0.5*dt/dx * (F[:,0,:] - F[:,1,:])

    return U_ebev 



def time_step_size(U, cfl, dx):
    """ 
        Function that calculates the size of the time step

        Args:
            U: ndarray(M+2,3); the conserved elements 2D array for the time step 
            cfl: float; courant number 
            dx: float; mesh size

        Returns:
            dt: float; the size of the time step 
    """
    c_ratio = 1.4
    
    # Converting the entire vector into the prim variables for ease of calculation
    W = np.empty(U.shape)
    for i, val in enumerate(U):
        W[i] = cmprin.consv_prim_1D(val)

    speeds = np.abs(W[:,1]) + np.sqrt(c_ratio * W[:,2] / W[:,0])

    dt = cfl * dx / np.max(speeds)

    return dt



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
    ebev_L, ebev_R = U_ebev[:,0], U_ebev[:,1] 

    # Forming the flux vector
    intercell_flux = np.empty(((ebev_R.shape[0]-1),3))
    for i in range(M+1):
        intercell_flux[i] = hllc.hllc_1D(ebev_R[i], ebev_L[i+1])

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
    time_evolved_U[1:-1] = U[1:-1] + dt/dx * (intercell_flux[:-1] - intercell_flux[1:])

    # Applying boundary conditions
    time_evolved_U = boundary_conditions(time_evolved_U)

    return time_evolved_U
    



def plotter(U, count, time):
    """ 
        Plots stuff

        Args:
            U: ndarray(M+2,3); the conserved elements 2D array for the time step 

    """
    
    # Obtaining vector in form of the primitive variables for ease of plotting
    W = np.empty(U.shape)
    for i, val in enumerate(U):
        W[i] = cmprin.consv_prim_1D(val)

    X = np.linspace(0,1,num=U.shape[0]-2)

    # Plotting variables
    plt.figure()
    plt.suptitle("Time t = " + str(time))

    plt.subplot(311)
    plt.plot(X, W[1:-1,0], 'bo-', linewidth=0.25, markersize=0.25)
    plt.title("Density vs Position")
    plt.ylabel("Density (kg/m^3)")
    plt.xlabel("Position (m)")

    plt.subplot(312)
    plt.plot(X, W[1:-1,1], 'bo-', linewidth=0.25, markersize=0.25)
    plt.title("Velocity vs Position")
    plt.ylabel("Velocity (m/s^2)")
    plt.xlabel("Position (m)")
    
    plt.subplot(313)
    plt.plot(X, W[1:-1,2], 'bo-', linewidth=0.25, markersize=0.25)
    plt.title("Pressure vs Position")
    plt.ylabel("Pressure (Pa)")
    plt.xlabel("Position (m)")

    plt.tight_layout()
    plt.savefig('test-##_itn-' + str(count) + '.png')
#    plt.show()

    # Plotting Residuals
    # TODO



def main():

    # Obtain DATA from the user
    print("\n\nEnter the following data")
    cfl = float(input("CFL Number: "))
    M = int(input("Number of data points: "))
    time = float(input("Time to which the solution should be obtained: "))
    print("\nThe length of the data line is assumed to be 1m. Starts from 0 and ends at 1m")
    plot_freq = int(input("\nFrequency of plotting (based on number of time steps): "))


    # Initialising some values 
    dx = 1/M

    # Creating the arrays given the user input 
    U = np.zeros((M+2, 3))

    # Initialising the array
    U = initialisation(U, M)

    # Applying the boundary conditions
    U = boundary_conditions(U)

    # Calling the iteration step (marching in time)
    count = 1
    march_time = 0

    while (march_time <= time):
        if ( count <= 5 ):
            dt = time_step_size(U, 0.2, dx)
        else: 
            dt = time_step_size(U, cfl, dx)
        march_time += dt

        # Obtaining conserved variable vector at march_time 
        U = stepping_time(U, dt, dx, M) 
        
        # Producing Plots
        if ( count % plot_freq == 0 ):
            plotter(U, count, march_time)
        




if __name__ == "__main__":
    main()
