using PyPlot

export plot_trajectories, plot_deviation, plot_trajectories_and_deviation,
		plot_responses, plot_responses_timeslice, plot_hatR, plot_hatR_diag, 
		plot_hatTheta

#font_size = 16
#rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
#font0 = Dict(
#        "font.size" => font_size,
#        "axes.labelweight" => "normal",
#        "axes.labelsize" => font_size,
#        "xtick.labelsize" => font_size,
#        "ytick.labelsize" => font_size,
#        "legend.fontsize" => font_size,
#		)
#merge!(rcParams, font0)


# Plot trajectories
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


# Plot Deviation
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


# Plot trajectories and deviation
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


# Plot responses
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


# Plot responses timeslice
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


# Plot hatR
function plot_hatR(tspan, hatR::Array{Float64,3}; quadR=false)

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
		if !quadR
			axes[i].set_title(L"$\hat{R}_1(\tau,\tau -\tau^\prime)$"*" for species $(i)")
		else 
			axes[i].set_title(L"$\hat{R}_2(\tau,\tau -\tau^\prime)$"*" for species $(i)")
		end
		i == ncols ? fig.colorbar(img, ax=vec(axes)) : nothing
	end
	
end


# Plot hatR diagional
function plot_hatR_diag(tspan, hatR::Array{Float64,3})

	num_species = size(hatR,1)
	dt = tspan[2] - tspan[1]
	
    figure(figsize=(20,5))
    for i in 1:num_species
        plot(tspan[1:end-1], [hatR[i,t+1,t]*dt for t in 1:length(tspan)-1])
    end
    legend(1:num_species)
    xlabel("Time")
    ylabel(L"$\hat{R}_1(\tau,\tau_-)\,\Delta t$")
    
end


# Plot hatTheta
function plot_hatTheta(plf::Plefka, fields::Fields)

	num_species = size(fields.hatTheta1,1)
	figure(figsize=(20,5))
	
	if plf.alphaOrder == 1
	
		for i in 1:num_species
			plot(fields.time[1:end-1], fields.hatTheta1[i,1:end-1])
		end
		legend(1:num_species)
		xlabel("Time")
		ylabel(L"$\tilde{\theta}_1$")
	
	elseif plf.alphaOrder == 2
		
		subplot(121)
		for i in 1:num_species
			plot(fields.time[1:end-1], fields.hatTheta1[i,1:end-1])
		end
		legend(1:num_species)
		xlabel("Time")
		ylabel(L"$\tilde{\theta}_1$")
		
		subplot(122)
		for i in 1:num_species
			plot(fields.time[1:end-1], fields.hatTheta2[i,1:end-1])
		end
		legend(1:num_species)
		xlabel("Time")
		ylabel(L"$\tilde{\theta}_2$")
	
	end
	
end
