using IterTools

export c_mn, e_i, m_list, n_list

"""
	c_mn(m, n, s_i, r_i, k3, y)
	
	Calculates the c_mn for general networks
	
	Uses as input
	m		Vector of integer
	n		Vector of integer
	s_i		stoichiometric coefficient for products
	r_i		stoichiometric coefficient for reactants
	k3		reaction rate for interactions
	y		molecule numbers
	
	Returns the time dependent c_mn  

"""
function c_mn(m::Any, n::Any, s_i::Matrix{Int64}, r_i::Matrix{Int64}, k3::Vector{Float64}, y::Matrix{Float64})

	num_species = size(y,1)		# number of species
	num_int = length(k3)		# number of interaction reactions
	len_ts = size(y,2)			# length of time grid
	res = zeros(len_ts)			# return vector
	prod_mu = ones(len_ts)		# product with mu_i
	prod1, prod2 = 1, 1			# initialize products
	
	for β in 1:num_int
		
		for i in 1:num_species
		
			prod1 *= binomial(s_i[β,i], m[i])
			prod2 *= binomial(r_i[β,i], m[i])
			
			prod_mu *= binomial(r_i[β,i], n[i])
			
			if r_i[β,i] > n[i]
				prod_mu .*=  y[i,:].^(r_i[β,i]-n[i])
			end
				
		end
		
		res += k3[β] * (prod1 - prod2) * prod_mu
		
	end
	
	return res

end


"""
	e_i(n, i)
	
	Returns a unit vector of length n with
	a 1 at position i

"""
function e_i(n::Int64, i::Int64)

	if n<=0 || i<=0 || i>n
		error("Invalid unit vector!")
	end
	
	e_i = zeros(Int64, n)
	e_i[i] = 1
	
	return e_i
end


"""
	mn_list(max_mn)
	
	Generic function to calculate either m_list() or n_list()

"""
function mn_list(max_mn::Vector{Int64})
	
	num_species = length(max_mn)
	
	x = [collect(0:max_mn[j]) for j in 1:num_species]
	m_list = collect(Iterators.product(x...))
	
	return reshape(m_list, length(m_list), 1)

end


"""
	n_list(r_i)

	Returns a list of possible vectors n
	The input r_i is the reactant coefficient matrix
	
"""
function n_list(r_i)
	max_n = maximum(r_i, dims=1)
	return mn_list(vec(max_n))
end


"""
	m_list(s_i, r_i)
	
	Returns a list of possible vectors m
	The input s_i and r_i are the stoichiometric 
	matrices
	
"""
function m_list(s_i, r_i)
	max_m = maximum([maximum(r_i, dims=1); maximum(s_i, dims=1)], dims=1)
	return mn_list(vec(max_m))
end

