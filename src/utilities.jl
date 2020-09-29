export Result, Responses, Fields, Plefka, status_update

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

mutable struct Fields_lin2
	hatTheta1::Array{Float64,2}
	hatTheta2::Array{Float64,2}
end

mutable struct Fields_quad1
	hatTheta1::Array{Float64,2}
	hatR1::Array{Float64,2}
end

mutable struct Fields_quad2
	hatTheta1::Array{Float64,2}
	hatTheta2::Array{Float64,2}
	hatR1::Array{Float64,2}
	hatR2::Array{Float64,3}
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
