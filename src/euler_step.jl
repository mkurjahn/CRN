using LinearAlgebra

export euler_step

#include("../examples/abc/abc_update_fields.jl")
#include("../examples/gene/gene_update_fields.jl")


"""
    euler_step(x0, params, tspan, plf::Plefka, s_i, r_i)
    
    Initializes the fields and runs the dynamics. 
    Uses an Euler step integrator.
    
    Returns a Struct with time and trajectories


"""
function euler_step(x0, params, tspan, plf::Plefka, s_i, r_i)

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
		fields = update_fields(plf, s_i, r_i, params[3], i, dt, y, resp, fields)
		
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
	update_fields(plf::Plefka, s_i, r_i, k3, t_i, dt, y, resp, fields)
	
	Updates the fields according to the theoretic equations 
	for general reaction networks
	
	Returns the updated fields

"""
function update_fields(plf::Plefka, s_i, r_i, k3, t_i, dt, y, resp, fields)
	
	N = size(y,1)		# number of species
	len_ts = size(y,2)	# length of time grid
	z = zeros(Int, N)
	
	for i in 1:N 
	
	if plf.orderParam == "linear"
			
		# update Theta1	
		fields.hatTheta1[i,t_i] = -c_mn(e_i(N,i), z, s_i, r_i, k3, y)[t_i]

		if plf.alphaOrder == 2
			
			# update Theta2
			sum_n = zeros(len_ts)
			for n in n_list(r_i)
				if sum(n) > 1
					prod_resp = ones(len_ts)
					for j in 1:N
						prod_resp .*= factorial(n[j])*(resp[j,t_i,:].^n[j])
					end
					sum_n += c_mn(e_i(N,i), n, s_i, r_i, k3, y)[t_i]*c_mn(n, z, s_i, r_i, k3, y).*prod_resp
				end
			end
			fields.hatTheta2[i,t_i] = -2*dt*sum(sum_n)
		
		end
	
	elseif plf.orderParam == "quad"
	
		# update Theta1 and hatR1
		c_eiei = c_mn(e_i(N,i), e_i(N,i), s_i, r_i, k3, y)[t_i]
		fields.hatTheta1[i,t_i] = -c_mn(e_i(N,i), z, s_i, r_i, k3, y)[t_i] + c_eiei*y[i,t_i]
		fields.hatR1[i,t_i+1,t_i] = -c_eiei/dt
		
		if plf.alphaOrder == 2
			
			# update Theta2 and hatR2
			sum_T2 = zeros(len_ts)
			sum_R2 = zeros(len_ts)
			local_R2 = 0
			for n in n_list(r_i)
				if sum(n) > 1
					prod_resp = ones(len_ts)
					for j in 1:N
						prod_resp .*= factorial(n[j])*(resp[j,t_i,:].^n[j])
					end
					c_n0 = c_mn(n, z, s_i, r_i, k3, y)
					c_nei = c_mn(n, e_i(N,i), s_i, r_i, k3, y)
					c_ein = c_mn(e_i(N,i), n, s_i, r_i, k3, y)[t_i]
					c_einei = c_mn(e_i(N,i), n .+ e_i(N,i), s_i, r_i, k3, y)[t_i]
					sum_T2 += (c_ein*c_n0 .- c_einei*c_n0*y[i,t_i] .- c_ein*c_nei.*y[i,:]).*prod_resp
					sum_R2 += c_ein*c_nei.*prod_resp
					local_R2 += sum(c_einei*c_n0.*prod_resp)
				elseif n[i] != 1
					prod_resp = ones(len_ts)
					for j in 1:N
						prod_resp .*= factorial(n[j])*(resp[j,t_i,:].^n[j])
					end
					c_nei = c_mn(n, e_i(N,i), s_i, r_i, k3, y)
					c_ein = c_mn(e_i(N,i), n, s_i, r_i, k3, y)[t_i]
					sum_T2 -= c_ein*c_nei.*y[i,:].*prod_resp
					sum_R2 += c_ein*c_nei.*prod_resp
				end
			end
			fields.hatTheta2[i,t_i] = -2*dt*sum(sum_T2)
			fields.hatR2[i,t_i+1,t_i] = -2*local_R2
			fields.hatR2[i,t_i,:] = -2*sum_R2
			
		end
		
	end
	
	end
	
	return fields

end


"""
	update_responses(plf::Plefka, resp, params, delta_t, t_i, hatR1, hatR2)
	
	Updates the responses
	
	Returns the updated response functions
	
"""
function update_responses(plf::Plefka, resp, params, delta_t, t_i, hatR1, hatR2)

	if plf.orderParam == "linear"
	
		resp[:,t_i+1,:] = (1 .- params[2]*delta_t).*resp[:,t_i,:]
	
	elseif plf.orderParam == "quad"
		
		if plf.alphaOrder == 1
			for k in t_i:-1:1
				x = plf.α*delta_t^2*sum((hatR1[:,t_i+1,:].*resp[:,:,k])[:,k:t_i], dims=2)
				resp[:,t_i+1,k] = (1 .- params[2]*delta_t).*resp[:,t_i,k] .- x
			end
			
		elseif plf.alphaOrder == 2
		
			for k in t_i:-1:1
				x = plf.α*delta_t^2*sum((hatR1[:,t_i+1,:].*resp[:,:,k])[:,k:t_i], dims=2)
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
	
		dydt = params[1] .- params[2] .* y[:,t_i] .- plf.α*fields.hatTheta1[:,t_i] .- plf.α*delta_t*sum(y[:,:].*fields.hatR1[:,t_i+1,:], dims=2)
		
		if plf.alphaOrder == 2
			dydt .-= 0.5*plf.α^2*fields.hatTheta2[:,t_i] .+ 0.5*plf.α^2*delta_t*sum(y[:,:].*fields.hatR2[:,t_i,:], dims=2)
		end
		
	end
	
	return dydt
end

	
