
def newton_raphson_method():
	start_step_value = float(input("Enter the value for the initial approximation: "))

	tolerance = float(input("Enter the error tolerance for the method: "))
	iteration_limit = float(input("Enter the iteration limit for the method: "))

	count = 0

	approximate_solution = newton_raphson_iteration_step(start_step_value, tolerance, count, iteration_limit)
	if ( approximate_solution != None ):
		print("The approximate solution for the zero of the given function by Newton Raphson Method is: ")
		print(approximate_solution)



def newton_raphson_iteration_step(step_value, tolerance, count, iteration_limit):
	updated_step_value = step_value - f(step_value) / f_dash(step_value)

	count += 1 
	if ( count > iteration_limit ):
		print("\nIteration limit has been exceeded. Solution unable to converge.")
		return None

	if ( (abs(updated_step_value - step_value) <= tolerance) | (f(updated_step_value) == 0)):
		return updated_step_value
	else: 
		return newton_raphson_iteration_step(updated_step_value, tolerance, count, iteration_limit)



def f(x):
	# User defines the function x here

def f_dash(x):
	# User defines the derivative of x here 



newton_raphson_method()
