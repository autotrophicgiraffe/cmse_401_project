# Level set function
function level(point)
	x = point[1]
	y = point[2]
	return x^2+y^2-1.0^2;
end

# greadient of the level set function


# dirichlet boundary conditions
function dirichlet_bc(point)
	return 0.0
end

# source function to make 0.5-exp(-20.0*((x-0.25)^2+(y-0.5)^2))-exp(-20.0*((x-0.75)^2+(y-0.5)^2))
# as a manufactured solution
function source_function(point,time)
	x = point[1]
	y = point[2]
	#return 8.0*(x^2+y^2)-1.0
	return exp(-time)*(0.75-7.0*(x^2+y^2))
end

function exact_solution(point,time)
    x = point[1]
    y = point[2]
    r = sqrt(x^2+y^2)
    nu = 7
    k77 = 31.4227941922
    return besselj(nu,k77*r)*cos(7.0*atan(y,x))*cos(k77*time)
end

function initial_condition(point)
    x = point[1]
    y = point[2]
    r = sqrt(x^2+y^2)
    nu = 7
    k77 = 31.4227941922
    return besselj(nu,k77*r)*cos(7.0*atan(y,x))
    #return exp(-100*(x^2+y^2))
end

function coefficientbeta(x,y)
	beta=1.0#0.25-(x^2+y^2)
	return beta
end
