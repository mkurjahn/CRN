using Pkg
cd(@__DIR__)
Pkg.activate("../../")
using CRN
using PyPlot
PyPlot.PyDict(PyPlot.matplotlib."rcParams")["font.size"] = 16
using Statistics

function params_study()

	# Parameters
	num_species = 3     # number of species
	num_int = 3         # number of interaction reaction

	# Reaction constants
	k1 = [8.0, 8.0, 3.0]      # Creation
	k2 = [3.0, 2.0, 1.5]    # Annihiliation
	k3 = [1.2, 0.3, 0.5]    # Interaction
	k = [k1, k2, k3]

	# Stoichiometric
	s_i = zeros(Int, num_int, num_species)
	r_i = copy(s_i)
	s_i[1,:] = [1 1 0]
	r_i[1,:] = [1 0 0]
	s_i[2,:] = [0 0 1]
	r_i[2,:] = [1 1 0]
	s_i[3,:] = [1 1 0]
	r_i[3,:] = [0 0 1]

	# Times
	t_init = 0.0        # Start time
	t_final = 5.0       # End time
	delta_t = 0.01      # time step
	tspan = collect(t_init:delta_t:t_final)

	# Initial condition
	x0 = [3., 1., 2.]
	
	# Plefka
	plf_sim_21 = Plefka(1.0, 2, "linear")
	plf_sim_22 = Plefka(1.0, 2, "quad")
	
	# Return
	r_gil = []
	r_pl1 = []
	r_pl2 = []
	
	# Decrease of k1,k2
	d = [1.0, 0.5, 0.1, 0.05, 0.01]
	
	for i in 1:length(d)
	
		println("Progress: run ", i)
		
		k1_new = k1 .* d[i]
		k2_new = k2 .* d[i]
		k2_new[2] = 1.0
		k = [k1_new, k2_new, k3]
		
		res_gil = gillespie_avg_v2(x0, k, s_i, r_i, tspan, 200000)
		res_plf_21 = euler_step(x0, k, tspan, plf_sim_21, s_i, r_i)
		res_plf_22 = euler_step(x0, k, tspan, plf_sim_22, s_i, r_i)
		
		push!(r_gil, res_gil)
		push!(r_pl1, res_plf_21[1])
		push!(r_pl2, res_plf_22[1])
	
	end
	
	return [r_gil, r_pl1, r_pl2]
	
end


function plot_study(r)

	(r_gil, r_pl1, r_pl2) = r
	
	for i in 1:length(r_gil)
	
		plot_deviation(r_pl1[i], r_gil[i])
		savefig("dev_plf1_gil_$(i).png")
		
		plot_deviation(r_pl2[i], r_gil[i])
		savefig("dev_plf2_gil_$(i).png")
		
		plot_deviation(r_pl1[i], r_pl2[i])
		savefig("dev_plf_plf_$(i).png")
	
	end

end


function steady_states(r)

	(r_gil, r_pl1, r_pl2) = r
	
	# num of results (gil+plf1+plf2) × num of params study × num of species
	S_states = zeros(length(r), length(r_gil), length(r_gil[1].data[:,1]))
	
	for i in 1:length(r)
		for j in 1:length(r_gil)
			S_states[i,j,:] = mean(r[i][j].data[:,end-10:end], dims=2)
		end
	end
	
	return S_states

end


function plot_steady_states(S_states)

	c = ["red", "blue", "green"]
	for i in 1:size(S_states, 3)	# num of species
		for j in 1:size(S_states, 1)	# num of results
			plot(1:size(S_states,2), S_states[j,:,i], color=c[j])
		end
		i == 1 ? legend(["GIL", "PLF_lin", "PLF_quad"]) : nothing
	end
end
