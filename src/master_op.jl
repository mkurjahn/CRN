using IterTools
using LinearAlgebra
using Distributions
using OrdinaryDiffEq

export get_max_num, master_operator, steady_state_masterOP, dynamics_masterOP,
		dynamics_masterOP_diag


"""
	get_max_num(x0::Vector{Float64}, precision::Float64)
	
	Calculates the vector of maximum number of species allowed where the
	Poisson distribution is below a certain probability (precision)
	
	Mean of Poisson distribution for each species is the initial condition x0
	
	Returns a vector of max_num

"""
function get_max_num(x0::Vector{Float64}, precision::Float64)
	
	num_species = length(x0)
	max_num = zeros(Int, num_species)
	for i in 1:num_species
		x,n = 1,0
		while abs(x) > precision
			x = 1-cdf(Poisson(x0[i]), n)
			n += 1
		end
		max_num[i] = n
	end
	return max_num
	
end
"""
	master_operator(p::Parameters, max_num; α=1.0)

	Calculates the master operator for the reaction network
	where the state space will be truncated at a maximum number 
	of species `max_num`
	
	Input
	p		Parameters struct for the reaction network
	max_num	number of species where the state space is truncated
	α		Plefka expansion parameter, default is α=1.0
	
	Returns the master operator and state space in a tuple
	
"""
function master_operator(p::Parameters, max_num::Vector{Int}; α=1.0)
	
	num_species = length(p.k[1])	# number of species
	params = p.k	# reaction rates
	s_i = p.s_i		# Stoichiometric products
	r_i = p.r_i		# Stoichiometric reactants
	
	# state space
	x = [collect(0:max_num[j]) for j in 1:num_species]
	state_space = reduce(vcat, collect(Iterators.product(x...)))
	
	# master operator
	len_state_space = prod(max_num .+ 1)
	master = zeros(len_state_space, len_state_space)
	
	for s in state_space
	
		# get index from s (assume only one appearence)
		idx_s = findfirst(x->x==s, state_space)
		
		for j in 1:num_species
			
			# Annihiliation
			t = collect(s)
			t[j] += 1
			t = tuple(t...)
			if t[j] <= max_num[j]
				idx_t = findfirst(x->x==t, state_space)		# index of t
				master[idx_s, idx_t] += params[2][j]*state_space[idx_t][j]
			end
			
			# Creation
			t = collect(s)
			t[j] -= 1
			t = tuple(t...)
			if t[j] >= 0
				idx_t = findfirst(x->x==t, state_space) 	# index of t
				master[idx_s, idx_t] += params[1][j]
			end
			
		end
			
		# Interaction
		for β in 1:length(params[3])
			
			new_state = tuple([(s[m] - s_i[β,m] + r_i[β,m]) for m in 1:num_species]...)
			
			if all(state_space[1] .<= new_state) && all(state_space[end] .>= new_state)
				idx = findfirst(x->x==new_state, state_space)
				r = α*params[3][β]
				for n in 1:num_species
					for p in 0:r_i[β,n]-1
						r *= (state_space[idx][n] - p)
					end
				end
				master[idx_s, idx] += r
			end
		
		end
	
	end
	
	# sum along columns must be zero
	master[diagind(master)] = -sum(master, dims=1)
	
	return (master, state_space)

end


function master_operator(p::Parameters, max_num::Int; α=1.0)
	max_num_vec = max_num .* ones(Int, length(p.k[1]))
	return master_operator(p::Parameters, max_num_vec; α=α)
end


"""
	steady_state_masterOP(master, state_space)
	
	Calculates the steady state of the master operator and state space
	
	Returns the steady state

"""
function steady_state_masterOP(master, state_space, E::Eigen)

	num_species = length(state_space[1])	# number of species
	
	# index of largest eigenvalue
	max_eig_idx = argmax(real.(E.values))
	
	# sum of eigenvectors of the highest eigenvalue
	sum_eigvecs = sum(abs.(E.vectors[:,max_eig_idx]))
	
	x0 = zeros(num_species)		# returned steady state
	
	for i in 1:length(state_space)
		for j in 1:num_species
			x0[j] += state_space[i][j] * abs(E.vectors[i,max_eig_idx]) / sum_eigvecs
		end
	end
	
	return x0

end


function steady_state_masterOP(master, state_space)
	E = eigen(master, sortby=nothing)		# eigenvalues and eigenvectors
	return steady_state_masterOP(master, state_space, E)
end

"""
	calc_mean_masterOP(p, state_space)
	
	Calculates the mean concentration with the probability vector `p`
	and state space
	
	Returns the calculated mean

"""
function calc_mean_masterOP(p::Vector{Float64}, state_space)

	num_species = length(state_space[1])	# number of species
	x = zeros(num_species)					# returned mean
	
	# sum of probabilities, should be == 1
	sum_p = sum(abs.(p))
	
	for i in 1:length(state_space)
		for j in 1:num_species
			x[j] += state_space[i][j] * abs(p[i]) / sum_p
		end
	end
	
	return x

end


"""
	initial_distr_state_space(state_space, x0)
	
	Calculates the initial Poisson distribution with means `x0` 
	for the different species.
	
	Returns the initial probability distribution vector

"""
function initial_distr_state_space(state_space, x0::Vector{Float64})

	num_species = length(state_space[1])	# number of species
	p0 = ones(length(state_space))			# initial prob. distr.
	
	for i in 1:length(state_space)
		for j in 1:num_species
			p0[i] *= pdf(Poisson(x0[j]), state_space[i][j])
		end
	end
	
	return p0./sum(p0)

end


"""
	dynamics_masterOP(master, state_space, tspan, x0, method)
	
	Runs the dynamics for the master operator using an euler step integrator.
	
	Input
	master		master operator
	state_space	truncated state space
	tspan		time grid
	x0			initial condition for copy numbers
	method		integration algorithm. If not specified: default is Tsit5()
	
	You can find a list of available solvers at https://diffeq.sciml.ai/stable/solvers/ode_solve/
	
	Returns a Result struct with time_grid and mean copy numbers

"""
function dynamics_masterOP(master, state_space, tspan::Vector{Float64}, x0::Vector{Float64}, method; rtn_prob=false)

	# Differential equation
	function f!(du,u,M,t)
		du .= M*u
	end
	
	# Solve the differential equation
	p0 = initial_distr_state_space(state_space, x0)    
    P = ode_solver(f!, master, tspan, p0, method; savedt=true)
    
    # Convert probability distributions to mean copy numbers
	y = zeros(length(x0), length(tspan))    
	for t in 1:length(tspan)
		y[:,t] = calc_mean_masterOP(P.data[:,t], state_space)
	end
	
	if rtn_prob
		return (P.data, Result(tspan, y))
	else
		return Result(tspan, y)
	end
end


function dynamics_masterOP(master, state_space, tspan::Vector{Float64}, x0::Vector{Float64})
	println("No integration algorithm specified! Using Tsit5() as default!")
	dynamics_masterOP(master, state_space, tspan, x0, Tsit5)
end


"""
	dynamics_masterOP_diag(master, state_space, E::Eigen, tspan, x0)
	
	Runs the dynamics for the master operator by diagonalizing the master
	operator and directly calculates the matrix exponential.
	
	Solves the master equation   dP/dt = M*P   by the solution
	
		P(t) = exp(M*t) * P0
	
	Input
	master		master operator
	state_space	truncated state space
	E			Eigen values and vectors, calculate by `eigen(master)`
	tspan		time grid
	x0			initial condition for copy numbers
	
	Returns a Result struct with time_grid and mean copy numbers

"""
function dynamics_masterOP_diag(master, state_space, E::Eigen, tspan, x0; rtn_prob=false)

	er = E.vectors
	el = inv(er)
	d = E.values #real.(E.values)
	
	p0 = initial_distr_state_space(state_space, x0) 
	y = zeros(length(x0), length(tspan))
	p = zeros(length(state_space), length(tspan))
	
	elp0 = el*p0
	
	for j in 1:length(tspan)
		for i in 1:length(state_space)
			@inbounds p[i,j] = abs(sum(er[i,:].*exp.(d*tspan[j]).*elp0))
		end
		p[:,j] ./= sum(p[:,j])
		y[:,j] = CRN.calc_mean_masterOP(p[:,j], state_space)
	end
	
	if rtn_prob
		return (p, Result(tspan, y))
	else
		return Result(tspan, y)
	end

end


function dynamics_masterOP_diag(master, state_space, tspan, x0)
	E = eigen(master, sortby=nothing)		# eigenvalues and eigenvectors
	return dynamics_masterOP_diag(master, state_space, E, tspan, x0)
end

