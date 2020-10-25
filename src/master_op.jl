using IterTools
using LinearAlgebra

#include("/home/maxi/ownCloud/Master/Simulations/testing/initialize.jl")

function master_operator(max_num, num_species, p::Parameters, α)

	params = p.k
	s_i = p.s_i
	r_i = p.r_i
	
	x = [collect(0:max_num) for j in 1:num_species]
	state_space = reduce(vcat, collect(Iterators.product(x...)))
	
	master = zeros((max_num+1)^num_species, (max_num+1)^num_species)
	
	for s in state_space
	
		# get index from s (assume only one appearence)
		idx_s = findfirst(x->x==s, state_space)
		
		for j in 1:num_species
			
			# Creation 
			t = collect(s)
			t[j] += 1
			t = tuple(t...)
			if t[j] <= max_num
				idx_t = findfirst(x->x==t, state_space) # index of t
				master[idx_s, idx_t] += params[2][j]*state_space[idx_t][j]
			end
			
			# Annihilation
			t = collect(s)
			t[j] -= 1
			t = tuple(t...)
			if t[j] >= 0
				idx_t = findfirst(x->x==t, state_space) # index of t
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
	
	master[diagind(master)] = -sum(master, dims=1)
	
	return (master, state_space)

end


function steady_state_master_op(master, state_space)

	num_species = length(state_space[1])
	
	E = eigen(master, sortby=nothing)
	
	max_eig_idx = argmax(real.(E.values))
	sum_eigvecs = sum(abs.(E.vectors[:,max_eig_idx]))
	#return eig_val
	x0 = zeros(num_species)
	
	for i in 1:length(state_space)
		for j in 1:num_species
			x0[j] += state_space[i][j] * abs(E.vectors[i,max_eig_idx]) / sum_eigvecs
		end
	end
	
	return x0

end


#p = ABC_params()
#(m,s) = master_operator(10, 3, p, 1.0)
#steady_state_master_op(m, s)

