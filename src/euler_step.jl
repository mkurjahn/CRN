using LinearAlgebra

export euler_step



"""
    euler_step(x0, params, tspan, plf::Plefka, s_i, r_i, invFunc)
    
    Initializes the fields and runs the dynamics. 
    Uses an Euler step integrator.
    
    Uses as input
    x0		Initial condition
    params	Reaction rates
    tspan	time grid vector
    plf		Plefka struct
    s_i,r_i Stoichiometric matrix, reactions in rows
    invFunc	Function, on which the expansion is based.
    
    If no function `invFunc` is specified, the default f(x)=1+x function
    will be used. Note that the defined function has to have the same
    Taylor expansion as f(x)=1+x.
    
    Returns a 3-element list with trajectories, responses and fields		


"""
function euler_step(x0, params, tspan, plf::Plefka, s_i, r_i, invFunc)

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
	
	# Initialize fields (depending on Plefka theory)
	fields = initialize_fields(plf, num_species, len_ts)	
	
	# run dynamics
	for i in 1:len_ts-1
		
		# update fields
		if plf.alphaOrder == 1
			update_fields!(fields, s_i, r_i, params[3], i, dt, y)
		else
			update_fields!(fields, s_i, r_i, params[3], i, dt, y, resp)
		end
		
		# update response functions
		if plf.orderParam == "linear"
			update_responses!(resp, params, dt, i)
		else
			update_responses!(resp, params, dt, i, plf.α, fields, invFunc)
		end
		
		# Euler step
		y[:,i+1] = y[:,i] + dt * y_derivative(y, params, i, dt, plf.α, fields, invFunc)
		
	end
	
	trajectories = Result(tspan, y)
	responses = Responses(tspan, resp)
	
	return [trajectories, responses, fields]
	
end

function euler_step(x0, params, tspan, plf::Plefka, s_i, r_i)
	@. f_theta(x) = 1 + x
	return euler_step(x0, params, tspan, plf, s_i, r_i, f_theta)
end

function euler_step(p::Parameters, plf::Plefka)
	return euler_step(p.x0, p.k, tspan(p), plf, p.s_i, p.r_i)
end

function euler_step(p::Parameters, plf::Plefka, f_theta)
	return euler_step(p.x0, p.k, tspan(p), plf, p.s_i, p.r_i, f_theta)
end


"""
	initialize_fields(plf::Plefka, num_species, len_ts)
	
	Initializes the fields depending on the Plefka theory
	
	Returns a struct with the needed fields

"""
function initialize_fields(plf::Plefka, num_species, len_ts)

	if plf.orderParam == "linear" && plf.alphaOrder == 1
		hatTheta1 = zeros(num_species, len_ts)
		fields = Fields_lin1(hatTheta1)

	elseif plf.orderParam == "linear" && plf.alphaOrder == 2
		hatTheta1 = zeros(num_species, len_ts)
		hatTheta2 = zeros(num_species, len_ts)
		fields = Fields_lin2(hatTheta1, hatTheta2)

	elseif plf.orderParam == "quad" && plf.alphaOrder == 1
		hatTheta1 = zeros(num_species, len_ts)
		hatR1 = zeros(num_species, len_ts)
		fields = Fields_quad1(hatTheta1, hatR1)
		
	elseif plf.orderParam == "quad" && plf.alphaOrder == 2
		hatTheta1 = zeros(num_species, len_ts)
		hatTheta2 = zeros(num_species, len_ts)
		hatR1 = zeros(num_species, len_ts)
		hatR2 = zeros(num_species, len_ts, len_ts)
		fields = Fields_quad2(hatTheta1, hatTheta2, hatR1, hatR2)
	end
	
	return fields

end


"""
	update_fields!(fields::Fields_lin1, s_i, r_i, k3, t_i, dt, y)
	
	Updates the fields according to the theoretic equations 
	for general reaction networks in the linear-1 case
	
	Returns the updated fields

"""
function update_fields!(fields::Fields_lin1, s_i, r_i, k3, t_i, dt, y)
	
	N = size(y,1)		# number of species
	len_ts = size(y,2)	# length of time grid
	z = zeros(Int, N)
	
	for i in 1:N 
		fields.hatTheta1[i,t_i] = -c_mn(e_i(N,i), z, s_i, r_i, k3, y, t_i)
	end

end


"""
	update_fields!(fields::Fields_lin2, s_i, r_i, k3, t_i, dt, y, resp)
	
	Updates the fields according to the theoretic equations 
	for general reaction networks in the linear-2 case
	
	Returns the updated fields

"""
function update_fields!(fields::Fields_lin2, s_i, r_i, k3, t_i, dt, y, resp)
	
	N = size(y,1)		# number of species
	len_ts = size(y,2)	# length of time grid
	z = zeros(Int, N)
	
	for i in 1:N 
				
		# update Theta1	
		fields.hatTheta1[i,t_i] = -c_mn(e_i(N,i), z, s_i, r_i, k3, y, t_i)

		# update Theta2
		sum_n = zeros(len_ts)
		for n in n_list(r_i)
			if sum(n) > 1
				prod_resp = calc_prod_resp(n, resp[:,t_i,:])
				sum_n += c_mn(e_i(N,i), n, s_i, r_i, k3, y, t_i)*c_mn(n, z, s_i, r_i, k3, y).*prod_resp
			end
		end
		fields.hatTheta2[i,t_i] = -2*dt*sum(sum_n)
		
	end

end


"""
	update_fields!(fields::Fields_quad1, s_i, r_i, k3, t_i, dt, y)
	
	Updates the fields according to the theoretic equations 
	for general reaction networks in the quad-1 case
	
	Returns the updated fields

"""
function update_fields!(fields::Fields_quad1, s_i, r_i, k3, t_i, dt, y)
	
	N = size(y,1)		# number of species
	len_ts = size(y,2)	# length of time grid
	z = zeros(Int, N)
	
	for i in 1:N 
		c_eiei = c_mn(e_i(N,i), e_i(N,i), s_i, r_i, k3, y, t_i)
		fields.hatTheta1[i,t_i] = -c_mn(e_i(N,i), z, s_i, r_i, k3, y, t_i) + c_eiei*y[i,t_i]
		fields.hatR1[i,t_i] = -c_eiei/dt	
	end

end


"""
	update_fields!(fields::Fields_quad2, s_i, r_i, k3, t_i, dt, y, resp)
	
	Updates the fields according to the theoretic equations 
	for general reaction networks in the quad-2 case
	
	Returns the updated fields

"""
function update_fields!(fields::Fields_quad2, s_i, r_i, k3, t_i, dt, y, resp)
	
	N = size(y,1)		# number of species
	len_ts = size(y,2)	# length of time grid
	z = zeros(Int, N)
	
	for i in 1:N 
	
		# update Theta1 and hatR1
		c_eiei = c_mn(e_i(N,i), e_i(N,i), s_i, r_i, k3, y, t_i)
		fields.hatTheta1[i,t_i] = -c_mn(e_i(N,i), z, s_i, r_i, k3, y, t_i) + c_eiei*y[i,t_i]
		fields.hatR1[i,t_i] = -c_eiei/dt

		# update Theta2 and hatR2
		sum_T2 = zeros(len_ts)
		sum_R2 = zeros(len_ts)
		local_R2 = 0
		for n in n_list(r_i)
			if sum(n) > 1
				prod_resp = calc_prod_resp(n, resp[:,t_i,:])
				c_n0 = c_mn(n, z, s_i, r_i, k3, y)
				c_nei = c_mn(n, e_i(N,i), s_i, r_i, k3, y)
				c_ein = c_mn(e_i(N,i), n, s_i, r_i, k3, y, t_i)
				c_einei = c_mn(e_i(N,i), n .+ e_i(N,i), s_i, r_i, k3, y, t_i)
				sum_T2 += (c_ein*c_n0 .- c_einei*c_n0*y[i,t_i] .- c_ein*c_nei.*y[i,:]).*prod_resp
				sum_R2 += c_ein*c_nei.*prod_resp
				local_R2 += sum((c_einei*c_n0.*prod_resp)[1:t_i])
			elseif sum(n) == 1 && n[i] != 1 
				prod_resp = calc_prod_resp(n, resp[:,t_i,:])
				c_nei = c_mn(n, e_i(N,i), s_i, r_i, k3, y)
				c_ein = c_mn(e_i(N,i), n, s_i, r_i, k3, y, t_i)
				sum_T2 -= c_ein*c_nei.*y[i,:].*prod_resp
				sum_R2 += c_ein*c_nei.*prod_resp
			end
		end
		fields.hatTheta2[i,t_i] = -2*dt*sum(sum_T2)
		fields.hatR2[i,t_i+1,t_i] = -2*local_R2
		fields.hatR2[i,t_i,:] = -2*sum_R2
		
	end
	
end


"""
	update_responses!(resp, params, delta_t, t_i)
	
	Updates the responses for linear order parameters
	
	Returns the updated response functions
	
"""
function update_responses!(resp, params, delta_t, t_i)

	resp[:,t_i+1,:] = (1 .- params[2]*delta_t).*resp[:,t_i,:]
	resp[:,t_i+1,t_i+1] .= 1.0
	
end


"""
	update_responses!(resp, params, delta_t, t_i, α, fields::Fields_quad1, invFunc)
	
	Updates the responses for quadratic order parameters (alpha lin)
	
	Returns the updated response functions
	
"""
function update_responses!(resp, params, delta_t, t_i, α, fields::Fields_quad1, invFunc)

	for i in 1:size(resp,1)
		resp[i,t_i+1,:] = (1 - params[2][i]*delta_t - α*delta_t^2*fields.hatR1[i,t_i]) .* resp[i,t_i,:]
	end
	resp[:,t_i+1,t_i+1] .= 1.0
	
end


"""
	update_responses!(resp, params, delta_t, t_i, α, fields::Fields_quad2, invFunc)
	
	Updates the responses for quadratic order parameters (alpha quad)
	
	Returns the updated response functions
	
"""
function update_responses!(resp, params, delta_t, t_i, α, fields::Fields_quad2, invFunc)

	for i in 1:size(resp,1)
		for k in t_i:-1:1
			x = α*delta_t^2*fields.hatR1[i,t_i]*resp[i,t_i,k]
			xx = 0.5*α^2*delta_t^2*sum((fields.hatR2[i,t_i,:].*resp[i,:,k])[k:t_i])
			x == 0 ? corr = 0 : corr = x * invFunc(xx/x)
			resp[i,t_i+1,k] = (1 - params[2][i]*delta_t)*resp[i,t_i,k] - corr
		end
	end
	resp[:,t_i+1,t_i+1] .= 1.0

end



"""
	y_derivative(y, params, t_i, delta_t, α, fields::Fields_lin1, invFunc)
	
	calculates the dy/dt derivative for the euler integration
	
	returns dydt
	
"""
function y_derivative(y, params, t_i, delta_t, α, fields::Fields_lin1, invFunc)
	
	dydt = params[1] .- params[2] .* y[:,t_i] .- α*fields.hatTheta1[:,t_i]
	return dydt
	
end


"""
	y_derivative(y, params, t_i, delta_t, α, fields::Fields_lin2, invFunc)
	
	calculates the dy/dt derivative for the euler integration
	
	returns dydt
	
"""
function y_derivative(y, params, t_i, delta_t, α, fields::Fields_lin2, invFunc)
	
	dydt = zeros(size(y,1))
	for i in 1:size(y,1)
		t1 = α*fields.hatTheta1[i,t_i]
		t2 = 0.5*α^2*fields.hatTheta2[i,t_i]
		t1 == 0.0 ? x = 0 : x = t1 * invFunc(t2/t1)
		dydt[i] = params[1][i] - params[2][i] * y[i,t_i] - x
	end
	return dydt
		
end



"""
	y_derivative(y, params, t_i, delta_t, α, fields::Fields_quad1, invFunc)
	
	calculates the dy/dt derivative for the euler integration
	
	returns dydt
	
"""
function y_derivative(y, params, t_i, delta_t, α, fields::Fields_quad1, invFunc)
	
	t1 = α*fields.hatTheta1[:,t_i]
	r1 = α*delta_t*y[:,t_i].*fields.hatR1[:,t_i]
	dydt = params[1] .- params[2] .* y[:,t_i] .- t1 .- r1
	return dydt
	
end



"""
	y_derivative(y, params, t_i, delta_t, α, fields::Fields_quad2, invFunc)
	
	calculates the dy/dt derivative for the euler integration
	
	returns dydt
	
"""
function y_derivative(y, params, t_i, delta_t, α, fields::Fields_quad2, invFunc)
	
	dydt = zeros(size(y,1))
	for i in 1:size(y,1)
		t1 = α*fields.hatTheta1[i,t_i]
		t2 = 0.5*α^2*fields.hatTheta2[i,t_i]
		r1 = α*delta_t*y[i,t_i].*fields.hatR1[i,t_i]
		r2 = 0.5*α^2*delta_t*sum(y[i,:].*fields.hatR2[i,t_i,:])
		t1+r1 == 0.0 ? x = 0 : x = (t1+r1) * invFunc((t2+r2)/(t1+r1))
		dydt[i] = params[1][i] - params[2][i] * y[i,t_i] - x
	end
	return dydt

end


	
