"""
    Python scripts that applies Godunov Method to solve a system of linear, hyperbolic
    equations

    System of Equations is of the form:
        U_t + A * U_x = 0
            
            where U_t is the partial derivatives of U wrt t 
                U_x is the partial derivative of U wrt x 
                A is the constant coefficient matrix (we shall make lots of references
                this is in documentation of the script)
"""
# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electrical & Electronics Engg. Minor 



import numpy as np 
import math as m 
import scipy.integrate as integrate



def initiali_tsar(length, num_data_pts, U):
    """
        Function that initialises the values for Godunov's Method at t = 0

        Args:
            length - float ; length of the data line 
            num_data_pts - int ; number of data points we want to sample for

        Return: 
            U - ndarray(number_of_eqautions, num_data_pts + 2) ; contains the 
                average value 
    """ 

    # ADD definitions of the initialisation functions for each element in U 
    
    # Creating a list of initialisation functions, each element of which is 
    # a call to a function 
    init_func_list = [] 

    dx = length/(num_data_pts + 1)
    
    # Initialising the values at t = 0
    # Note that we are using the finite volume method
    for init, i in enumerate(init_func_list):
        for j in range(1, num_data_pts + 1):
            U[i,j] = (1/dx) * integrate.quad(init, j-0.5, j+0.5)
    
    U[:,0], U[:,-1] = U[:,1], U[:,-2]

    return U
    
    
def splitting_A(A):
    """ 
        Function that takes the coefficient matrix A and obtains 
        the diagonal eigenvalue matrix, Lambda; the matrix of right
        eigenvectors, K; the Lambda_plus matrix; the Lambda_minus 
        matrix; 
        
        Args: 
            A - 2D ndarray ; constant coefficient matrix of the sys of 
                linear hyperbolic equations 
        
        Returns:
            A_plus - 2D ndarray ; (ref. text)
            A_minus - 2D ndarray ; (ref. text)
            max_characterisitic - float ; the value of the largest eigenvalue
                of the constant coefficient matrix A 
    """

    eigval, eigvec = np.linalg.eig(A)
    max_characterisitic = np.max(eigval)
    
    Lambda_plus = np.zeros(eigval.shape[0])
    Lambda_minus = np.zeros(eigval.shape[0])
    for eigvals, i in enumerate(eigval):
        Lambda_plus[i] = max(eigvals, 0)
        Lambda_minus[i] = min(eigvals, 0)

    A_plus =  np.dot(eigvec, np.dot(np.diag(Lambda_plus), np.linalg.inv(eigvec)))
    A_minus =  np.dot(eigvec, np.dot(np.diag(Lambda_minus), np.linalg.inv(eigvec)))

    return A_plus, A_minus, max_characterisitic


def godunov_method_sys_linear_hyperbolic():
    """
        This is the main function that uses Godunov's method to obtain a solution 
        of the system of linear, hyperbolic equations over the time interval 
        specified by the user. 
        The user also provides the length of the region to obtain solutions for and 
        the number of data points (the precision to which) he wants the solution to be 
        obtained for
        This function also produces plots after a defined number of time steps
    """

    length = float(input("\nEnter the length of the data line: ")) 
    num_data_pts = int(input("Enter the number of data (the precision to which) the 
                                solution should be obtained for: "))
    time_period = float(input("\nEnter the time period for which the solution should
                                be obtained"))
    cfl = float(input("CFL coefficient: "))
    print("\nUser should note that the value incremental value for each time step will 
                                be obtained by applying the CFL stability condition."))
        
    # EDIT ACCORDINGLY 
    # ADD defintion of the constant coefficient matrix A 
    A = np.array()
    
    # Creating the matrix U 
    U = np.zeros((A.shape[0], num_data_pts + 2))
    
    # Initialising the values of U 
    U = initiali_tsar(length, num_data_pts, U)
    
   
    # GODUNOV's METHOD 
    A_plus, A_minus, max_characterisitic = splitting_A(A) 
    flux = np.zeros(num_data_pts + 1)
    
    dx = lenth / (num_data_pts + 1)
    t = 0
    dt = cfl * dx / max_characterisitic
    
    while ( t <= time_period ):
        for flx, i in enumerate(flux):
            flx = np.dot(A_plus, U[:,i]) + np.dot(A_minus, U[:,i+1])
        
        # The Godunov Update Step
        U[:, 1:-1] += (dt/dx) * ( flux[:,:-1] - flux[:,1:] )
        # Applying the boundary condition 
        U[:, 0], U[:, -1] = U[:, 1], U[:, -2] 
        
        # ADD plotting feature  
