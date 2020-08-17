using OrdinaryDiffEq

export ode_solver


"""
    ode_solver(f!, k, tspan)
    
    Solves the differential equations for the mass actions
    and Plefka problem. Vector k contains the reactions constants 
    and tspan is a tupel with start and end time.
    
    Returns the result in a struct with time and mean copy 
    number of molecules. 


"""
function ode_solver(f!, k, tspan, dt, x0)
	
	prob = ODEProblem(f!, x0, tspan, k)

	sol = solve(prob, Tsit5(), saveat=dt)

	col = []
	for i in 1:length(x0)
	    push!(col, map(A->A[i], sol.u))
	end
    
    xar = reshape(cat(col..., dims=1), length(sol.t), length(x0))
    return Result(sol.t, xar')
end


function ode_solver(f!, k, tspan, x0)
	return ode_solver(f!, k, tspan, 0.01, x0)
end
