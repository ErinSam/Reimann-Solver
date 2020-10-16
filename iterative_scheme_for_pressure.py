# Python Implementation of Iterative Scheme for Finding the Pressure 
# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electrical & Electronics Minor 

# Script uses the Newton-Raphson iterative procedure to find the root of the 
# pressure function, of the Reimann Problem with ideal gas Equation of State 

# The value of pressure in the star region can then be used to solve for the value 
# of velocity in the star region. From the velocity and pressure, we can then find 
# density and remaining unknowns using standard gas dynamics relations

#######################################################################################################
####################################################################################################### 


import numpy as np
import math as m
from single_variable_root_finding import *


def star_region_properties(**kwargs):
    
    # Ignore the poor implementation
    p_L = kwargs['p_L'] 
    p_R = kwargs['p_R'] 
    rho_L = kwargs['rho_L'] 
    rho_R = kwargs['rho_R'] 
    u_L = kwargs['u_L'] 
    u_R = kwargs['u_R'] 
    a_L = kwargs['a_L'] 
    a_R = kwargs['a_R'] 
    A_L = kwargs['A_L'] 
    A_R = kwargs['A_R'] 
    B_L = kwargs['B_L'] 
    B_R = kwargs['B_R'] 
    c_ratio = kwargs['c_ratio'] 


    # 3 Kinds of Initial Approximations
    # Two-Rarefaction Approximation 
    p_TR = pow( (a_L + a_R - 0.5*(c_ratio-1)*(u_R-u_L))/( (a_L/pow(p_L, (c_ratio-1)/(2*c_ratio))) +(a_R/pow(p_R, (c_ratio-1)/(2*c_ratio))) ), (2*c_ratio)/(c_ratio-1) )
    
    ## Approximation from a Linearised Solution based on Primitive Variables
    # p_PV = 0.5 * (p_L + p_R) - 0.125 * (u_R - u_L)*(rho_L + rho_R)*(a_L + a_R)
    # p_0 = max(TOL, p_PV) 
    ## Two-Shock Approximation 
    # g_R = m.sqrt(A_R/(p_0+B_R))
    # g_L = m.sqrt(A_L/(p_0+B_L))
    # p_TS = ( g_L*p_L + g_R*p_R - (u_R - u_L) )/( g_L + g_R )
    # p_00 = max(TOL, p_TS)
   

 
    # Calculating the pressure in the star region
    # Using Two-Rarefaction Approximation:
    p_star = newton_raphson_method(pressure_func, d_pressure_func, p_TR, p_L=p_L, p_R=p_R, rho_L=rho_L, rho_R=rho_R, u_L=u_L, u_R=u_R, a_L=a_L, a_R=a_R, A_L=A_L, A_R=A_R, B_L=B_L, B_R=B_R, c_ratio=c_ratio, calc="pressure")
    
    ## Using Approximation from a Linearised Solution based on Primitive Variables  
    # p_star_2 = nrm(pressure_func, d_pressure_func, p_0, p_L=p_L, p_R=p_R, rho_L=rho_L, rho_R=rho_R, u_L=u_L, u_R=u_R, a_L=a_L, a_R=a_R, A_L=A_L, A_R=A_R, B_L=B_L, B_R=B_R, c_ratio=c_ratio)
    ## Using Two-Shock Approximation
    # p_star_3 = nrm(pressure_func, d_pressure_func, p_00, p_L=p_L, p_R=p_R, rho_L=rho_L, rho_R=rho_R, u_L=u_L, u_R=u_R, a_L=a_L, a_R=a_R, A_L=A_L, A_R=A_R, B_L=B_L, B_R=B_R, c_ratio=c_ratio)

    # Calculating the velocity in the Star Region 
    u_star = 0.5 * (u_L + u_R) + 0.5 * pressure_func(p_star, p_L=p_L, p_R=p_R, rho_L=rho_L, rho_R=rho_R, u_L=u_L, u_R=u_R, a_L=a_L, a_R=a_R, A_L=A_L, A_R=A_R, B_L=B_L, B_R=B_R, c_ratio=c_ratio, calc="velocity") 
    
    # Calculating the densities in the Star Region
    # rho_star_L = rho_L * ( ( (p_star/p_L) + (c_ratio-1)/(c_ratio+1) )/( (c_ratio-1)/(c_ratio+1)*(p_star/p_L) + 1) )
    # rho_star_R = rho_R * ( ( (p_star/p_R) + (c_ratio-1)/(c_ratio+1) )/( (c_ratio-1)/(c_ratio+1)*(p_star/p_R) + 1) )

    
    # Printing the values
    # print("\n\nThe value of pressure in the star region, p*: ", p_star)
    # print("The value of velocity in the star region, u*: ", u_star)
    # print("The value of density in the left star region, rho*L: ", rho_star_L)
    # print("The value of density in the right star region, rho*R: ", rho_star_R)

    return p_star, u_star


# Defining the pressure function 
def pressure_func(p, **kwargs):

    # Ignore the poor implementation
    p_L = kwargs['p_L'] 
    p_R = kwargs['p_R'] 
    rho_L = kwargs['rho_L'] 
    rho_R = kwargs['rho_R'] 
    u_L = kwargs['u_L'] 
    u_R = kwargs['u_R'] 
    a_L = kwargs['a_L'] 
    a_R = kwargs['a_R'] 
    A_L = kwargs['A_L'] 
    A_R = kwargs['A_R'] 
    B_L = kwargs['B_L'] 
    B_R = kwargs['B_R'] 
    c_ratio = kwargs['c_ratio'] 
    
    # function L is given by 
    def f_L():
        if ( p > p_L ):
            return ( (p - p_L) * m.sqrt( A_L/(p+B_L) ) )
        else:
            return ( (2*a_L)/(c_ratio-1) * ( pow(p/p_L, (c_ratio-1)/(2*c_ratio)) - 1 ) )
    
    # function R is given by 
    def f_R():
        if ( p > p_R ):
            return ( (p - p_R) * m.sqrt( A_R/(p+B_R) ) )
        else:
            return ( (2*a_R)/(c_ratio-1) * ( pow(p/p_R, (c_ratio-1)/(2*c_ratio)) - 1) )

    if ( kwargs['calc'] == "pressure" ):    
        return ( f_L() + f_R() + (u_R - u_L) )
    elif ( kwargs['calc'] == "velocity" ):
        return ( f_R() - f_L() )



# Defining first dervative of the pressure function
def d_pressure_func(p, **kwargs):

    # Ignore the poor implementation
    p_L = kwargs['p_L'] 
    p_R = kwargs['p_R'] 
    rho_L = kwargs['rho_L'] 
    rho_R = kwargs['rho_R'] 
    u_L = kwargs['u_L'] 
    u_R = kwargs['u_R'] 
    a_L = kwargs['a_L'] 
    a_R = kwargs['a_R'] 
    A_L = kwargs['A_L'] 
    A_R = kwargs['A_R'] 
    B_L = kwargs['B_L'] 
    B_R = kwargs['B_R'] 
    c_ratio = kwargs['c_ratio'] 

    # first derivative of function L is given by 
    def d_f_L():
        if ( p > p_L ):
            return ( m.sqrt(A_L/(p+B_L)) * (1 - (p - p_L)/(2*(B_L + p))) )
        else:
            return ( 1/(rho_L*a_L) * pow(p/p_L, -(c_ratio+1)/(2*c_ratio)) )

    # first derivative of function R is given by 
    def d_f_R():
        if ( p > p_R ):
            return ( m.sqrt(A_R/(p+B_R)) * (1 - (p - p_R)/(2*(B_R + p ))) )
        else:
            return ( 1/(rho_R*a_R) * pow(p/p_R, -(c_ratio+1)/(2*c_ratio)) )

    return ( d_f_L() + d_f_R() )  



#######################################################################################################
#######################################################################################################

# star_region_properties()

