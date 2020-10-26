using OrdinaryDiffEq

export ode_solver


"""
    ode_solver(f!, k, ts::Vector{Float64}, x0::Vector{Float64}, method; savedt=true)
    
    Solves differential equations!
    
    Input
    f!		set of differential equations
    k		reactions rates for the network
    ts		time grid vector
    x0		initial condition
    method	integration algorithm. If not specified: default is Tsit5()
	
	You can find a list of available solvers at https://diffeq.sciml.ai/stable/solvers/ode_solve/
	
	savedt	If true, the result will be saved after each time step.
			If false, an adaptive time step is used.
			For method==Euler, this has to be true and dt will be set to the
			discretization in the time grid vector `ts`
    
    Returns the result in a struct with time and mean copy 
    numbers of molecules. 


"""
function ode_solver(f!, k, ts::Vector{Float64}, x0::Vector{Float64}, method; savedt=true)
	
	tspan = (ts[1], ts[end])
	prob = ODEProblem(f!, x0, tspan, k)
	
	if method==Euler
		println("Euler algoritm needs choice of dt, setting savedt=true")
		delta_t = (ts[2]-ts[1])
		sol = solve(prob, Euler(), dt=delta_t, saveat=delta_t)
	elseif savedt
		sol = solve(prob, method(), saveat=(ts[2]-ts[1]))
	else
		sol = solve(prob, method())
	end
	
	col = []
	for i in 1:length(x0)
	    push!(col, map(A->A[i], sol.u))
	end
    
    xar = reshape(cat(col..., dims=1), length(sol.t), length(x0))
    return Result(sol.t, xar')
end


function ode_solver(f!, k, tspan::Tuple, x0)
	println("ode_solver: This is a old version! Please use a vector for tspan
	and a specified method for the integration algorithm (default: Tsit5()). ")
	ts = [tspan[1], tspan[end]]
	ode_solver(f!, k, ts, x0, Tsit5; savedt=false)
end

function ode_solver(f!, k, ts::Vector{Float64}, x0::Vector{Float64})
	println("No integration algorithm specified! Using Tsit5() as default!")
	ode_solver(f!, k, ts, x0, Tsit5)
end
