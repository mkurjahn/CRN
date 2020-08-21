using LinearAlgebra

export euler_step

include("../examples/abc/abc_update_fields.jl")


"""
    euler_step(x0, params, tspan, plf::Plefka)
    
    Initializes the fields and runs the dynamics. 
    Uses an Euler step integrator.
    
    Returns a Struct with time and trajectories


"""
function euler_step(x0, params, tspan, plf::Plefka)

	global num_species = length(params[1])
	
	len_ts = length(tspan)
	dt = tspan[2] - tspan[1]
	
	# Initialize trajectory
	y = zeros(num_species, len_ts)
	y[:,1] = x0
	
	# Initialize response functions
	resp = zeros(num_species, len_ts, len_ts)
	for i in 1:num_species
		resp[i,:,:] = Matrix{Int}(I, len_ts, len_ts)
	end
	
	# Initialize fields (even if we don't need them)
	hatTheta1 = zeros(num_species, len_ts)
	hatTheta2 = copy(hatTheta1)
	hatR1 = zeros(num_species, len_ts, len_ts)
	hatR2 = copy(hatR1)
	fields = Fields(tspan, hatTheta1, hatTheta2, hatR1, hatR2)
	
	
	# run dynamics
	for i in 1:len_ts-1
		
		# update fields
		fields = ABC_update_fields(plf, params[3][1], i, dt, y, resp, fields)
		
		# update response functions
		resp = update_responses(plf, resp, params, dt, i, fields.hatR1, fields.hatR2)
		
		# Euler step
		y[:,i+1] = y[:,i] + dt * y_derivative(plf, y, params, i, dt, fields)
		
	end
	
	trajectories = Result(tspan, y)
	responses = Responses(tspan, resp)
	
	return [trajectories, responses, fields]
	
end



"""
	update_responses(plf::Plefka, resp, params, delta_t, t_i, hatR1, hatR2)
	
	updates the responses hatR1 and hatR2
	
	returns the updated response functions
	
"""
function update_responses(plf::Plefka, resp, params, delta_t, t_i, hatR1, hatR2)

	if plf.orderParam == "linear"
	
		resp[:,t_i+1,:] = (1 .- params[2]*delta_t).*resp[:,t_i,:]
	
	elseif plf.orderParam == "quad"
		
		if plf.alphaOrder == 1
			for k in t_i:-1:1
				x = plf.α*delta_t^2*sum((hatR1[:,t_i,:].*resp[:,:,k])[:,k:t_i], dims=2)
				resp[:,t_i+1,k] = (1 .- params[2]*delta_t).*resp[:,t_i,k] .- x
			end
			
		elseif plf.alphaOrder == 2
		
			for k in t_i:-1:1
				x = plf.α*delta_t^2*sum((hatR1[:,t_i,:].*resp[:,:,k])[:,k:t_i], dims=2)
				xx = 0.5*plf.α^2*delta_t^2*sum((hatR2[:,t_i,:].*resp[:,:,k])[:,k:t_i], dims=2)
				resp[:,t_i+1,k] = (1 .- params[2]*delta_t).*resp[:,t_i,k] .- x .- xx
			end
		
		end
	
	end
	
	resp[:,t_i+1,t_i+1] .= 1.0	
	return resp
	
end



"""
	y_derivative(plf::Plefka, y, params, t_i, delta_t, hatTheta1, hatTheta2, hatR1, hatR2)
	
	calculates the dy/dt derivative for the euler integration
	
	returns dydt
	
"""
function y_derivative(plf::Plefka, y, params, t_i, delta_t, fields::Fields)
	
	dydt = zeros(num_species)
	
	if plf.orderParam == "linear"
	
		dydt = params[1] .- params[2] .* y[:,t_i] .- plf.α*fields.hatTheta1[:,t_i]

		if plf.alphaOrder == 2
			dydt .-= 0.5*plf.α^2*fields.hatTheta2[:,t_i]
		end
		
	elseif plf.orderParam == "quad"
	
		dydt = params[1] .- params[2] .* y[:,t_i] .- plf.α*fields.hatTheta1[:,t_i] .-
		plf.α*delta_t*sum(y[:,:].*fields.hatR1[:,t_i,:], dims=2)
		
		if plf.alphaOrder == 2
			dydt .-= 0.5*plf.α^2*fields.hatTheta2[:,t_i] .+
			0.5*plf.α^2*delta_t*sum(y[:,:].*fields.hatR2[:,t_i,:], dims=2)
		end
		
	end
	
	return dydt
end

	
