# Python Implementation of Newton-Raphon Method (copy with slight variation, to suit problem at hand, from my MA207 repo) 
# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electrical & Electronics Minor

########################################################################################################
########################################################################################################

def newton_raphson_method(f, f_dash, initial_approximation, **kwargs):

    # start_step_value = float(input("Enter the value for the initial approximation: "))
    # tolerance = float(input("\nEnter the error tolerance for root finding method: "))
    # iteration_limit = float(input("Enter the iteration limit for the method: "))
    start_step_value = initial_approximation
    tolerance = 1e-6
    iteration_limit = 20

    count = 0
    approximate_solution = newton_raphson_iteration_step(start_step_value, tolerance, count, iteration_limit, f, f_dash, **kwargs)
    if ( approximate_solution != None ):
        return approximate_solution


def newton_raphson_iteration_step(step_value, tolerance, count, iteration_limit, f, f_dash, **kwargs):
    updated_step_value = step_value - f(step_value, **kwargs) / f_dash(step_value, **kwargs)

    count += 1 
    if ( count > iteration_limit ):
        print("\nIteration limit has been exceeded. Solution unable to converge.")
        return None

    if ( (abs(updated_step_value - step_value) <= tolerance) | (f(updated_step_value, **kwargs) == 0)):
        return updated_step_value
    else: 
        return newton_raphson_iteration_step(updated_step_value, tolerance, count, iteration_limit, f, f_dash, **kwargs)

########################################################################################################
########################################################################################################
