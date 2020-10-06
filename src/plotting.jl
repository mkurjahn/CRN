using PyPlot

export plot_trajectories, plot_deviation, plot_trajectories_and_deviation,
		plot_responses, plot_responses_timeslice, plot_hatR, plot_hatR_diag, 
		plot_hatTheta, cut_res_tf, cut_resp_tf

"""
	cut_res_tf(res::Result, tf)
	
	Only time grid up to a specified final time tf which has to be
	a part of the vector res.time
	
	Returns the Result struct up to time tf

"""
function cut_res_tf(res::Result, tf)
	ts = findall(x->x==tf, res.time)
	length(ts) != 1 ? error("No matching time slice!") : t = ts[1]
	return Result(res.time[1:t], res.data[:,1:t])
end


"""
	cut_resp_tf(resp::Responses, tf)
	
	Only time grid up to a specified final time tf which has to be
	a part of the vector res.time
	
	Returns the Response struct up to time tf

"""
function cut_resp_tf(resp::Responses, tf)
	ts = findall(x->x==tf, resp.time)
	length(ts) != 1 ? error("No matching time slice!") : t = ts[1]
	return Responses(resp.time[1:t], resp.data[:,1:t,1:t])
end


"""
	plot_trajectories(res::Result)
	
	Plots the trajectories for all species
	
	Use plot_trajectories(res::Result, tf) for plotting up to final time tf

"""
function plot_trajectories(res::Result)

	num_species = size(res.data,1)
	figure(figsize=(20,5))
	
	for i in 1:num_species
		plot(res.time, res.data[i,:], linestyle="dashed")
	end
	
	legend(1:num_species)
	xlabel("Time")
	ylabel("Copy Numbers")
	
end

function plot_trajectories(res::Result, tf)
	plot_trajectories(cut_res_tf(res, tf))
end


"""
	plot_deviation(res::Result, ref::Result)
	
	Plots the deviation of res for all species to the reference ref 
	(in most cases e.g. Gillespie)
	
	deviation = | res - ref | / ref
	
	Use plot_deviation(res::Result, ref::Result, tf) for plotting up to
	final time tf

"""
function plot_deviation(res::Result, ref::Result)

	num_species = size(res.data,1)
	figure(figsize=(20,5))
	
	for i in 1:num_species
		plot(res.time, abs.(res.data[i,:] .- ref.data[i,:]) ./ ref.data[i,:], linestyle="dashed")
	end
	
	legend(1:num_species)
	xlabel("Time")
	ylabel("Relative deviation from reference")
	
end

function plot_deviation(res::Result, ref::Result, tf)
	plot_deviation(cut_res_tf(res, tf), cut_res_tf(ref, tf))
end


"""
	plot_trajectories_and_deviation(res::Result, ref::Result)
	
	Plots the trajectories and deviaton for all species.
	
	For more details, look at `plot_trajectories(...)`
	and `plot_deviation(...)`
	
	Use plot_trajectories_and_deviation(res::Result, ref::Result, tf) for 
	plotting up to final time tf

"""
function plot_trajectories_and_deviation(res::Result, ref::Result)

	num_species = size(res.data,1)
	figure(figsize=(20,5))
	
	subplot(121)
	for i in 1:num_species
		plot(res.time, res.data[i,:], linestyle="dashed")
	end
	legend(1:num_species)
	xlabel("Time")
	ylabel("Copy Numbers")
	
	subplot(122)
	for i in 1:num_species
		plot(res.time, abs.(res.data[i,:] .- ref.data[i,:]) ./ ref.data[i,:], linestyle="dashed")
	end
	legend(1:num_species)
	xlabel("Time")
	ylabel("Relative deviation from reference")
	
end

function plot_trajectories_and_deviation(res::Result, ref::Result, tf)
	plot_trajectories_and_deviation(cut_res_tf(res, tf), cut_res_tf(ref, tf))
end


"""
	plot_responses(resp::Responses)
	
	Plots the response functions for all species
	
	Use plot_responses(resp::Responses, tf) for plotting up to final time tf

"""
function plot_responses(resp::Responses)

	num_species = size(resp.data,1)
	
	ncols = 3
	vmin = minimum(resp.data)
	vmax = maximum(resp.data)
	lims = [resp.time[1], resp.time[end], resp.time[end], resp.time[1]]
	
	rr = copy(resp.data)

	for i in 1:num_species
		for j in 1:length(resp.time)
			rr[i,j,:] = circshift(resp.data[i,j,end:-1:1],j)
		end
	end

	fig, axes = subplots(figsize=(20,5), ncols=ncols, nrows=ceil(Int, num_species/ncols), sharey=true)
	
	for i in 1:num_species
		img = axes[i].imshow(rr[i,:,:], aspect="equal", extent=lims, vmin=vmin, vmax=vmax)
		axes[i].set_xlabel(L"$\tau - \tau^\prime$")
		axes[i].set_title(L"$ R(\tau,\tau -\tau^\prime) $"*" for species $(i)")
		i == ncols ? fig.colorbar(img, ax=vec(axes)) : nothing
	end
	
end

function plot_responses(resp::Responses, tf)
	plot_responses(cut_resp_tf(resp, tf))
end


"""
	plot_responses_timeslice(resp::Responses, timeslice)
	
	Plots the response for all species at a given timeslice
	
	R(t=timeslice, t')

"""
function plot_responses_timeslice(resp::Responses, timeslice)

	num_species = size(resp.data,1)
	ts = findall(x->x==timeslice, resp.time)
	length(ts) != 1 ? error("No matching time slice!") : t = ts[1]
	figure(figsize=(20,5))
	
	for i in 1:num_species
		plot(resp.time, resp.data[i,t,:], linestyle="dashed")
	end
	
	legend(1:num_species)
	xlabel("Time")
	ylabel(L"$ R(\tau="*string(timeslice)*L",\tau^\prime) $")

end


"""
	plot_hatR(tspan, fields::Fields_quad2)
	
	Plots the field hatR2 for the quadratic order parameter and alpha squared
	case. 
	
	Needs improvement if number of species is not divisible by 3.

"""
function plot_hatR(tspan, fields::Fields_quad2)

	hatR = fields.hatR2
	num_species = size(hatR,1)
	
	ncols = 3
	vmin = minimum(hatR)
	vmax = maximum(hatR)
	vmin == 0 && vmax == 0 ? vmax = 1.0 : nothing
	lims = [tspan[1], tspan[end], tspan[end], tspan[1]]
	
	rr = copy(hatR)

	for i in 1:num_species
		for j in 1:length(tspan)
			rr[i,j,:] = circshift(hatR[i,j,end:-1:1],j)
		end
	end

	fig, axes = subplots(figsize=(20,5), ncols=ncols, nrows=ceil(Int, num_species/ncols), sharey=true)
	
	for i in 1:num_species
		img = axes[i].imshow(rr[i,:,:], aspect="equal", extent=lims, vmin=vmin, vmax=vmax)
		axes[i].set_xlabel(L"$\tau - \tau^\prime$")
		axes[i].set_title(L"$\hat{R}_2(\tau,\tau -\tau^\prime)$"*" for species $(i)")
		i == ncols ? fig.colorbar(img, ax=vec(axes)) : nothing
	end
	
end


"""
	plot_hatR_diag(tspan, fields)
	
	Plots the hatR1 diagonal elements R1(τ,τ-Δt)

"""
function plot_hatR_diag(tspan, fields)

	hatR = fields.hatR1
	num_species = size(hatR,1)
	dt = tspan[2] - tspan[1]
	
    figure(figsize=(20,5))
    for i in 1:num_species
        plot(tspan[1:end-1], hatR[i,1:end-1]*dt)
    end
    legend(1:num_species)
    xlabel("Time")
    ylabel(L"$\hat{R}_1(\tau,\tau_-)\,\Delta t$")
    
end


# Plot hatTheta1
function plot_hatTheta1(tspan, fields)

	num_species = size(fields.hatTheta1,1)
	figure(figsize=(20,5))
		
	for i in 1:num_species
		plot(tspan[1:end-1], fields.hatTheta1[i,1:end-1])
	end
	legend(1:num_species)
	xlabel("Time")
	ylabel(L"$\tilde{\theta}_1$")
		
end

"""
	plot_hatTheta(tspan, fields::Fields_lin1)
	
	Plots the hat Theta field for linear order parameters and linear alpha

"""
function plot_hatTheta(tspan, fields::Fields_lin1)
	plot_hatTheta1(tspan, fields)
end

"""
	plot_hatTheta(tspan, fields::Fields_quad1)
	
	Plots the hat Theta field for quadratic order parameters and linear alpha

"""
function plot_hatTheta(tspan, fields::Fields_quad1)
	plot_hatTheta1(tspan, fields)
end


# Plot hatTheta1 and hatTheta2
function plot_hatTheta12(tspan, fields)

	num_species = size(fields.hatTheta1,1)
	figure(figsize=(20,5))
		
	subplot(121)
	for i in 1:num_species
		plot(tspan[1:end-1], fields.hatTheta1[i,1:end-1])
	end
	legend(1:num_species)
	xlabel("Time")
	ylabel(L"$\tilde{\theta}_1$")
	
	subplot(122)
	for i in 1:num_species
		plot(tspan[1:end-1], fields.hatTheta2[i,1:end-1])
	end
	legend(1:num_species)
	xlabel("Time")
	ylabel(L"$\tilde{\theta}_2$")
	
end

"""
	plot_hatTheta(tspan, fields::Fields_lin2)
	
	Plots the hat Theta fields 1 and 2 for linear order parameters and
	alpha squared

"""
function plot_hatTheta(tspan, fields::Fields_lin2)
	plot_hatTheta12(tspan, fields)
end

"""
	plot_hatTheta(tspan, fields::Fields_quad2)
	
	Plots the hat Theta fields 1 and 2 for quadratic order parameters and
	alpha squared

"""
function plot_hatTheta(tspan, fields::Fields_quad2)
	plot_hatTheta12(tspan, fields)
end
