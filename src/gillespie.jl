using Distributions
using Roots
using DataFrames
using StatsBase
using LinearAlgebra

export gillespie, gillespie_tspan, gillespie_full, gillespie_avg, gillespie_avg_v2


"""
    pfsample(w)
    
    Select a single random integer betwenn 1 and length(w) with 
    probabilities proportional to the weights given in w.
    
    Based on the function StatsBase.sample(w)
    
    
"""
function pfsample(w::Array{Float64,1})
    t = rand() * sum(w)
    n = length(w)
    i = 1
    cw = w[1]
    while cw < t && i < n
        i += 1
        @inbounds cw += w[i]
    end
    return i
end


"""
	calc_nu(num_species, s_i, r_i)
	
	Calculates the stochiometric matrix nu
	
	Uses as input
	num_species		number of species
	s_i				stochiometric matrix for product molecules
					(only interaction reaction)
	r_i				stochiometric matrix for reactant molecules 
					(only interaction reaction)
	
	Returns the stochiometric matrix for the systems inclusive the baseline
	creation and annihilation reaction
	

"""
function calc_nu(num_species::Int64, s_i::Matrix{Int64}, r_i::Matrix{Int64})
	eye = Matrix{Int}(I, num_species, num_species)
	return [eye; -eye; (s_i-r_i)]
end



"""
	prop_fct(x, k, r_i)
	
	Calculates the propensity function (reaction probability vector)
	
	Uses as input
	x		vector of molecule numbers
	params	reaction rate parameters
	r_i		stochiometric matrix for reactant molecules 
			(only interaction reaction)
	
	Returns the reaction probability vector

"""
function prop_fct(x, params::Vector{Vector{Float64}}, r_i::Matrix{Int64})
	num_int = length(params[3])
	Hint = ones(num_int)
	for m in 1:num_int
		for k in 1:length(params[1])
			for r in 0:r_i[m,k]-1
				@inbounds Hint[m] *= (x[k]-r)
			end
		end
	end
	return [params[1]; params[2].*x'; params[3].*Hint]
end



"""
	regular_timeGrid(res, ta)
	
	Takes the output from the irregular time discretization
	of the Gillespie algorithm and makes the time grid regular
	
	Uses as input
	res		Result struct from the gillespie function
	ta		Desired time discretization 
	
	Returns the evenly spaced time discretized Result
	
"""
function regular_timeGrid(res::Result, ta::Vector{Float64})
	
	num_species = size(res.data, 2)		# number of species
	l_ta = length(ta)					# length of time grid
	
	erg = zeros(l_ta, num_species)		# result 
	i = 1
	j = 1
	
	# loop
	while i <= l_ta
		if ta[i] <= res.time[j+1]
			erg[i,:] = res.data[j,:]
		else
			i -= 1
			j += 1
		end
		i += 1
	end

	return Result(ta, erg)
end



"""
    gillespie(x0, params, nu, r_i, tf)
    
    Gillespie simulation of one specific trajectory and
    irregular time discretization.
    
    Uses as input
    x0      Vector of molecule numbers
    params  Vector of reaction constants
    nu		Stoichiometric matrix, reactions in rows
    r_i		Reactant coefficients
    tf      Final time when the reaction has to stop
    
    Returns the result in a struct with time and copy 
    number of molecules. 


"""
function gillespie(x0::Vector{Int64}, params::Vector{Vector{Float64}}, nu::Matrix{Int64}, r_i::Matrix{Int64}, tf::Float64)
    
    ta = Vector{Float64}()
    t = 0.0
    push!(ta,t)

    # Set up initial x
    nstates = length(x0)
    x = copy(x0')
    xa = copy(x0)

    # Number of reactions
    numpf = size(nu, 1)

    nsteps = 0
    while t <= tf
        pf = prop_fct(x,params,r_i)
        sumpf = sum(pf)
        if sumpf == 0.0
            # break condition
            t += 2*tf
        end
        dt = rand(Exponential(1/sumpf))
        t += dt
        push!(ta,t)
        # update event
        ev = pfsample(pf)
        deltax = view(nu,ev,:)
        for i in 1:nstates
            @inbounds x[1,i] += deltax[i]
        end
        for xx in x
            push!(xa,xx)
        end
        # update nsteps
        nsteps += 1
    end
    xar = transpose(reshape(xa,length(x),nsteps+1))
    return Result(ta,xar)
end



"""
    gillespie_tspan(x0, params, nu, r_i, time_grid)
    
    Gillespie simulation of one specific trajectory and
    even spaced time discretization.
    
    Uses as input
    x0      Vector of molecule numbers
    params  Vector of reaction constants
    nu      Stoichiometric matrix, reactions in rows
    r_i		Reactant coefficients
    time_grid	Vector of time discretization
    
    Returns the result in a struct with time and copy 
    number of molecules. 


"""
function gillespie_tspan(x0::Vector{Int64}, params::Vector{Vector{Float64}}, nu::Matrix{Int64}, r_i::Matrix{Int64}, time_grid::Vector{Float64})

    # Set up initial x
    nstates = length(x0)
    x = copy(x0')
    xnew = copy(x)
    xa = Vector{Int64}()

    # Number of reactions
    numpf = size(nu, 1)

    tm = time_grid[1]
    tf = time_grid[end]
    
    for t in time_grid
    
		while tm <= t
		
			x = copy(xnew)
		    pf = prop_fct(x,params,r_i)
		    sumpf = sum(pf)
		    
		    if sumpf == 0.0
		    	# break condition
		        tm += 2*length(time_grid)
		    end
		    
		    # update event    
	        ev = pfsample(pf)
	        deltax = view(nu,ev,:)
	        for i in 1:nstates
	            @inbounds xnew[1,i] += deltax[i]
	        end

		    tm += rand(Exponential(1/sumpf))
			
		end
		
		for xx in x
			push!(xa,xx)
		end
	    
	end
		
	return transpose(reshape(xa,length(x),length(time_grid)))
	
end



"""
    gillespie_full(w0, params, s_i, r_i, ta, N)
    
    Full N Gillespie Simulations with saved copy numbers.
    
    Uses as input
    w0      Initial conditions drawn from a Poisson
            distribution with means w0
    params  Vector of reaction constants
    s_i,r_i Stoichiometric matrix, reactions in rows
    ta		Vector of time discretization
    num_runs   Number of runs/repetitions
    
    Returns the full Gillespie Simulation of N runs
     

"""
function gillespie_full(w0::Vector{Float64}, params::Vector{Vector{Float64}}, s_i::Matrix{Int64}, r_i::Matrix{Int64}, ta::Vector{Float64}, num_runs::Int64)

	# Stochiometric matrix
	nu = calc_nu(length(params[1]), s_i, r_i)
	
    # Array to store results
    copyN = zeros(Int32, length(ta), size(s_i,2), num_runs)
    
    # loop over runs
    @simd for i in 1:num_runs
	    x0 = rand.(Poisson.(w0))
	    erg = gillespie_tspan(x0, params, nu, r_i, ta)
	    copyN[:,:,i] = erg[:,:]
    end

	return copyN
end



"""
    gillespie_avg(w0, params, s_i, r_i, ta, N)
    
    Average of N Gillespie simulations.
    
    Same as gillespie_full, but averages over the N simulations.
    
    Returns the mean result in a struct with time and 
    mean copy number of molecules. 


"""
function gillespie_avg(w0::Vector{Float64}, params::Vector{Vector{Float64}}, s_i::Matrix{Int64}, r_i::Matrix{Int64}, ta::Vector{Float64}, num_runs::Int64)

    copyN = gillespie_full(w0, params, s_i, r_i, ta, num_runs)
    return Result(ta, reshape(mean(copyN, dims=3), length(ta), size(s_i,2))')

end



"""
	gillespie_avg_v2(w0, params, s_i, r_i, ta, num_runs)
	
	Average of N Gillespie simulations.
	
	Same as gillespie_avg, but the regular time grid is taken after
	the gillespie simulation and not within.
	
	Returns the mean result in a struct with time and
	mean copy number of molecules.

"""
function gillespie_avg_v2(w0::Vector{Float64}, params::Vector{Vector{Float64}}, s_i::Matrix{Int64}, r_i::Matrix{Int64}, ta::Vector{Float64}, num_runs::Int64)

	# Stochiometric matrix
	nu = calc_nu(length(params[1]), s_i, r_i)
	
    # Array to store results
    copyN = zeros(Int32, length(ta), size(s_i,2), num_runs)
    
    # loop over runs
    @simd for i in 1:num_runs
	    x0 = rand.(Poisson.(w0))
	    gil = gillespie(x0, params, nu, r_i, ta[end])
	    erg = regular_timeGrid(gil, ta).data
	    copyN[:,:,i] = erg[:,:]
    end

	return Result(ta, reshape(mean(copyN, dims=3), length(ta), size(s_i,2))')
	
end






