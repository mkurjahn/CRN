
export print_fields, eq_cmn, print_cmn, print_cmn_br


"""
	print_prod_resp(n, species)
	
	Prints the product of responses for
	n		Vector of integers (the n vector)
	species	Vector of strings with species names, e.g. ['A', 'B', 'C']
	
	No return value

"""
function print_prod_resp(n, species)

	prod_resp = ""
	for j in 1:length(n)
		if n[j] == 1
			prod_resp *= "·R"*species[j]
		elseif n[j] == 2
			prod_resp *= "·R"*species[j]*"²"
		elseif n[j] > 2
			prod_resp *= string(factorial(n[j]))
			prod_resp *= "·R"*species[j]*"^"*string(n[j])
		end
	end
	print(prod_resp)
	
end


"""
	eq_cmn(m, n, s_i::Matrix{Int64}, r_i::Matrix{Int64}, k, species; print=false)
	
	Calculates the binomial prefactor and product of means for the c_mn
	
	Input
	m		Vector of integers which are the 'm' values
	n		Vector of integers which are the 'n' values
	s_i,r_i	Stoichiometric matrices
	k		Vector of strings with reaction rates, e.g. ['k1', 'k2', 'k3']
	species	Vector of strings with species names, e.g. ['A', 'B', 'C']
	
	If print==false returns
	res		Vector of integers which is the prefactor of the c_mn
			and includes all binomial factors
	mus		Vector of strings with the products of means of the 
			different species according to the binomial factors
	
	Otherwise
	Call the function `print_cmn(res, mus)`
	
"""
function eq_cmn(m, n, s_i::Matrix{Int64}, r_i::Matrix{Int64}, k, species; print=false)
	
	num_int = length(k)
	res = zeros(Int, num_int)
	mus = ["" for i in 1:num_int]
		
	for β in 1:num_int
	
		prod_mu = 1					# product with mu_i
		prod1, prod2 = 1, 1			# initialize products
		
		mus[β] *= k[β] #*"·"
	
		for i in 1:length(species)
		
			prod1 *= binomial(s_i[β,i], m[i])
			prod2 *= binomial(r_i[β,i], m[i])
			prod_mu *= binomial(r_i[β,i], n[i])
			
			if r_i[β,i] == n[i]+1
				mus[β] *= "·μ"*species[i]
			elseif r_i[β,i] == n[i]+2
				mus[β] *=  "·μ"*species[i]*"²"
			elseif r_i[β,i] > n[i]+2
				mus[β] *=  "·μ"*species[i]*"^"*string(r_i[β,i]-n[i])
			end
			
		end
		
		res[β] += (prod1 - prod2) * prod_mu
	
	end
	
	if print
		print_cmn(res, mus)
	else
		return (res, mus)
	end
	
end


"""
	print_cmn(res::Vector{Int}, mus::Vector{String})
	
	Prints the c_mn
	
	Input
	res		Vector of integers which is the prefactor of the c_mn
			and includes all binomial factors
	mus		Vector of strings with the products of means of the 
			different species according to the binomial factors
	
	Usually takes as input the return values of `eq_cmn(...)`
	
	No return value

"""
function print_cmn(res::Vector{Int}, mus::Vector{String})
	
	if all(res .== 0)
		print("0")
	end
	
	firstprint = true
	for β in 1:length(res)
		if firstprint
			if res[β]==1 
				print(mus[β])
				firstprint=false
			elseif res[β]==-1
				print("-", mus[β])
				firstprint=false
			elseif abs(res[β]) > 1
				print(res[β], mus[β])
				firstprint=false
			end
		elseif abs(res[β]) == 1
			if res[β]>0 
				print("+ ", mus[β])
			elseif res[β]<0
				print("– ", mus[β])
			end
		elseif abs(res[β]) > 1
			if res[β]>0
				print("+ ", res[β], mus[β])
			elseif res[β]<0
				print("– ", -res[β], mus[β])
			end
		end
		if β != length(res) && res[β]*res[β+1] != 0
			print(" ")
		end
	end
	
end


"""
	print_cmn_br(res::Vector{Int}, mus::Vector{String})
	
	Same as print_cmn(res, mus) but with brackets before and after
	the expression

"""
function print_cmn_br(res::Vector{Int}, mus::Vector{String})
	print("(")
	print_cmn(res, mus)
	print(")")
end


"""
	print_fields(p::Parameters, plf::Plefka, k_param, species; disconnected_resp=true)
	
	Prints the fields according to the theoretical equations for linear 
	and/or quadratic order parameters and expansion order alpha and
	alpha^2, depending on the definition of plf::Plefka
	
	Input:
	p		Parameters for the reaction network
	plf		Plefka definition
	k_param	Vector of strings with reaction rates, e.g. ['k1', 'k2', 'k3']
	species	Vector of strings with species names, e.g. ['A', 'B', 'C']
	
	diconnected_resp	Boolean specifying if we use disconnected or connected 
						response functions for the calculations of the fields
	
	No output, only printing to the standard output  (stdout)


"""
function print_fields(p::Parameters, plf::Plefka, k_param, species; disconnected_resp=true)
	
	N = length(p.k[1])		# number of species
	z = zeros(Int, N)		# Vector of zeros
	c_mn2(m,n) = eq_cmn(m, n, p.s_i, p.r_i, k_param, species; print=false)
	
	if plf.orderParam == "linear"
	
		for i in 1:N
		
			println("––––––––––\nSpecies: ", i)
			
			# Theta1
			print("Theta1 = ")
			(r1,m1) = c_mn2(e_i(N,i), z)
			all(r1 .== 0) ? print("0") : print_cmn_br(-r1, m1)
			println()
				
			if plf.alphaOrder == 2
				
				firstprint = true
				print("Theta2 = ")
				for n in n_list(p.r_i)
					if sum(n) > 1
						(r3,m3) = c_mn2(e_i(N,i), n)
						(r4,m4) = c_mn2(n, z)
						if any(r3 .!= 0) && any(r4 .!= 0)
							firstprint ? print("-2Δt sum[ ") : print(" + ")
							print_cmn_br(r3, m3)
							print_cmn_br(r4, m4)
							print_prod_resp(n,species)
							firstprint = false
						end	
					end
				end
				firstprint ? println("0") : println(" ]")				
			end
			println()
		end
		
	elseif plf.orderParam == "quad"
	
		for i in 1:N
		
			println("––––––––––\nSpecies: ", i)
			
			# Theta1
			print("Theta1 = ")
			(r1,m1) = c_mn2(e_i(N,i), z)
			(r2,m2) = c_mn2(e_i(N,i), e_i(N,i))
			
			all(r1 .== 0) ? print("") : print_cmn_br(-r1, m1)
			if all(r1 .== 0) && all(r2 .== 0) 
				print("0")
			end
			if disconnected_resp
				print(" + ")
				print_cmn_br(r2, m2)
				println("·μ", species[i])
			else
				println()
			end
			
			# hatR1
			print("hatR1  = ")
			if all(r2 .== 0)
				println("0")
			else
				print_cmn_br(-r2, m2)
				println("/Δt")
			end
			
			if plf.alphaOrder == 2
			
				# Theta2
				firstprint = true
				print("Theta2 = ")
				
				for n in n_list(p.r_i)
				
					(r3,m3) = c_mn2(e_i(N,i), n)
					(r4,m4) = c_mn2(n, z)
					(r5,m5) = c_mn2(n, e_i(N,i))
					(r6,m6) = c_mn2(e_i(N,i), n .+ e_i(N,i))
						
					if sum(n) > 1
						
						all([r3; r4; r5; r6] .== 0) ? print("") : print(" + ")
						
						if any(r3 .!= 0) && any(r4 .!= 0)
							firstprint ? print("-2Δt sum[ ") : print(" + ")
							print_cmn_br(r3, m3)
							print_cmn_br(r4, m4)
							print_prod_resp(n,species)
							firstprint = false
						end
						if any(r6 .!= 0) && any(r4 .!= 0) && disconnected_resp
							firstprint ? print("-2Δt sum[ ") : print(" + ")
							print_cmn_br(-r6, m6)
							print_cmn_br(r4, m4)
							print("·μ", species[i])
							print_prod_resp(n,species)
						end
						if any(r3 .!= 0) && any(r5 .!= 0) && disconnected_resp
							firstprint ? print("-2Δt sum[ ") : print(" + ")
							print_cmn_br(-r3, m3)
							print_cmn_br(r5, m5)
							print("·μ", species[i])
							print_prod_resp(n,species)
						end
						if any([r3; r4; r5; r6] .!= 0)
							firstprint = false
						end
						
					
					elseif sum(n) == 1 && n[i] != 1
						
						if any(r3 .!= 0) && any(r5 .!= 0) && disconnected_resp
							firstprint ? print("-2Δt sum[ ") : print(" + ")
							print_cmn_br(-r3, m3)
							print_cmn_br(r5, m5)
							print("·μ", species[i])
							print_prod_resp(n,species)
							firstprint = false
						end
					
					end
					
				end
				
				firstprint ? println("0") : println(" ]")
				
				# hatR2
				firstprint = true
				print("hatR2  = ")
				
				for n in n_list(p.r_i)
				
					(r3,m3) = c_mn2(e_i(N,i), n)
					(r4,m4) = c_mn2(n, z)
					(r5,m5) = c_mn2(n, e_i(N,i))
					(r6,m6) = c_mn2(e_i(N,i), n .+ e_i(N,i))
						
					if sum(n) > 1
						
						if any(r6 .!= 0) && any(r4 .!= 0)
							firstprint ? print("-2·[ ") : print(" + ")
							print("sum( ")
							print_cmn_br(r6, m6)
							print_cmn_br(r4, m4)
							print("·δ(t,t'-)")
							print_prod_resp(n,species)
							print(" )")
							firstprint = false
						end
						if any(r3 .!= 0) && any(r5 .!= 0)
							firstprint ? print("-2·[ ") : print(" + ")
							print_cmn_br(r3, m3)
							print_cmn_br(r5, m5)
							print_prod_resp(n,species)
							firstprint = false
						end
						
					
					elseif sum(n) == 1 && n[i] != 1
						
						if any(r3 .!= 0) && any(r5 .!= 0)
							firstprint ? print("-2·[ ") : print(" + ")
							print_cmn_br(r3, m3)
							print_cmn_br(r5, m5)
							print_prod_resp(n,species)
							firstprint = false
						end
					
					end
					
				end
				
				firstprint ? println("0") :	println(" ]")
			
			end
			
			println()
			
		end
		
	end

end


