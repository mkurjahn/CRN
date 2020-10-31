import Base.display

export Result, Responses, Plefka, Parameters, tspan, status_update

"""
    struct Result
    
    Struct with the result of the chemical species.
    
    Contains:
    time    Vector of time steps
    data    Matrix of mean copy numbers for each species

"""
struct Result
    time::Vector{Float64}
    data::Array{Float64,2}
end

function Base.display(res::Result)
	display([res.time'; res.data])
end


"""
    struct Responses
    
    Struct with the result of the response functions.
    
    Contains:
    time    Vector of time steps
    data    3-dim Array with species, time, time

"""
struct Responses
    time::Vector{Float64}
    data::Array{Float64,3}
end

function Base.display(resp::Responses)
	display(permutedims(resp.data, (2,3,1)))
end


"""
	struct Fields
	
	Struct with all auxiliary fields 
	Differentiate between Plefka theory!
	
	Contains:
	hatTheta1, 
	hatTheta2, 
	hatR1, 
	hatR2

"""
mutable struct Fields_lin1
	hatTheta1::Array{Float64,2}
end
function Base.display(fields::Fields_lin1)
	display(fields.hatTheta1)
end

mutable struct Fields_lin2
	hatTheta1::Array{Float64,2}
	hatTheta2::Array{Float64,2}
end
function Base.display(fields::Fields_lin2)
	display(fields.hatTheta1)
	display(fields.hatTheta2)
end

mutable struct Fields_quad1
	hatTheta1::Array{Float64,2}
	hatR1::Array{Float64,2}
end
function Base.display(fields::Fields_quad1)
	display(fields.hatTheta1)
	display(fields.hatR1)
end

mutable struct Fields_quad2
	hatTheta1::Array{Float64,2}
	hatTheta2::Array{Float64,2}
	hatR1::Array{Float64,2}
	hatR2::Array{Float64,3}
end
function Base.display(fields::Fields_quad2)
	display(fields.hatTheta1)
	display(fields.hatTheta2)
	display(fields.hatR1)
	display(permutedims(fields.hatR2, (2,3,1)))
end


"""
	struct Plefka
	
	Struct with all information about the Plefka approximation.
	
	Contains:
	alpha		number between 0 and 1, typically == 1
	alphaOrder	order of alpha expansion, either 1 or 2
	orderParam	either 'linear' or 'quad'
	
"""
struct Plefka
	α::Float64
	alphaOrder::Int64
	orderParam::String
	function Plefka(α::Float64, alphaOrder::Int64, orderParam::String)
		if α < 0.0 || α > 1.0
			error("PLEFKA: invalid alpha value! Between 0 and 1!")
		end
		if alphaOrder != 1 && alphaOrder != 2
			error("PLEFKA: invalid alpha order! Either 1 or 2!")
		end
		if orderParam != "linear" && orderParam != "quad"
			error("PLEFKA: invalid order parameter! Either 'linear' or 'quad'!")
		end
		new(α, alphaOrder, orderParam)
	end
end


"""
	mutable struct Parameters
	
	inherits all important information about the 
	parameters of the system
	
	x0		initial condition
	k		reaction rates
	s_i		stoichiometric products
	r_i		stoichiometric reactants
	t_init	initial time
	t_final	final time
	delta_t	time discretization

"""
mutable struct Parameters
	x0::Vector{Float64}
	k::Vector{Vector{Float64}}
	s_i::Matrix{Int64}
	r_i::Matrix{Int64}
	t_init::Float64
	t_final::Float64
	delta_t::Float64
end


"""
	tspan(p::Parameters)
	
	return an array of the time grid
	according to the time discretization
	in the Parameters set p.

"""
function tspan(p::Parameters)
	return collect(p.t_init:p.delta_t:p.t_final)	
end


"""
	status_update(t, len_ts)
	
	Prints the status in percent from [0..100%] in 10% steps
	
	t		loop index
	len_ts	length of the loop
	
	Put in a loop!

"""
function status_update(t, len_ts)
	
	if t == 1
		println("[  0%]")
	elseif t%(len_ts÷10) == 0 && t < 0.91*len_ts
		println("[ ", round(Int, 100*t/len_ts), "%]")
	elseif t == len_ts-1
		println("[100%]")
	end

end
