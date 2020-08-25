"""

    gene_update_fields(plf::Plefka, k3, t_i, y, resp, hatTheta1, hatTheta2, hatR1, hatR2)
    
    Updates the fields for the gene example network. 
    Need to figure out how to generalize this!
    
    Input:
    plf		Plefka Struct
    k3		reaction rates for gene network
    t_i		time step
    y		trajectory for G_on, G_off, Protein
    resp	response functions
    hatTheta1, hatTheta2, hatR1, hatR2 are the auxiliary fields

	returns the updated fields (hatTheta1, hatTheta2, hatR1, hatR2)

"""
function gene_update_fields(plf::Plefka, k3, t_i, delta_t, y, resp, fields::Fields)

	if plf.orderParam == "linear"
		
		fields.hatTheta1[1,t_i] = k3[2]*y[1,t_i]*y[2,t_i] - k3[3]*y[3,t_i]
		fields.hatTheta1[2,t_i] = -k3[1]*y[1,t_i] + k3[2]*y[1,t_i]*y[2,t_i] - k3[3]*y[3,t_i]
		fields.hatTheta1[3,t_i] = -k3[2]*y[1,t_i]*y[2,t_i] + k3[3]*y[3,t_i]
	
		if plf.alphaOrder == 2
			
			x = 2*delta_t*k3[2]*sum( (k3[1]*y[1,:] .- k3[2]*y[1,:].*y[2,:] .+ k3[3]*y[3,:]) .* resp[1,t_i,:] .* resp[2,t_i,:] )
			fields.hatTheta2[1,t_i] = x
			fields.hatTheta2[2,t_i] = x
			fields.hatTheta2[3,t_i] = -x
		
		end
	
	elseif plf.orderParam == "quad"
	
		fields.hatTheta1[1,t_i] = -k3[3]*y[3,t_i]
		fields.hatTheta1[2,t_i] = -k3[1]*y[1,t_i] - k3[3]*y[3,t_i]
		fields.hatTheta1[3,t_i] = -k3[2]*y[1,t_i]*y[2,t_i]
		
		fields.hatR1[1,t_i+1,t_i] = k3[2]*y[2,t_i]/delta_t
		fields.hatR1[2,t_i+1,t_i] = k3[2]*y[1,t_i]/delta_t
		fields.hatR1[3,t_i+1,t_i] = k3[3]/delta_t
		
		if plf.alphaOrder == 2
			
			fields.hatTheta2[1,t_i] = -2*delta_t*k3[2]*( -sum(k3[3]*y[3,:].*resp[1,t_i,:].*resp[2,t_i,:]) + y[1,t_i]*sum((k3[1]*y[1,:] .- k3[2]*y[1,:].*y[2,:] .+ k3[3]*y[3,:]).*resp[2,t_i,:]) )
			fields.hatTheta2[2,t_i] = -2*delta_t*k3[2]*( -sum( (k3[1]*y[1,:] .+ k3[3]*y[3,:]).*resp[1,t_i,:].*resp[2,t_i,:] ) + y[2,t_i]*sum((k3[1]*y[1,:] .- k3[2]*y[1,:].*y[2,:] .+ k3[3]*y[3,:]).*resp[1,t_i,:]) )
			fields.hatTheta2[3,t_i] = -2*delta_t*k3[2]*sum( (k3[1]*y[1,:] .- k3[2]*y[1,:].*y[2,:]) .* resp[1,t_i,:] .* resp[2,t_i,:] )
			
			fields.hatR2[1,t_i,:] = 2*k3[2]*((k3[1] .- k3[2]*y[2,:]).*(resp[1,t_i,:] .+ y[1,:]) .+ k3[3]*y[3,:]).*resp[2,t_i,:]
			fields.hatR2[2,t_i,:] = 2*k3[2]*(-k3[2]*y[1,:].*(resp[2,t_i,:] .+ y[2,:]) .+ k3[1]*y[1,:] .+ k3[3]*y[3,:]).*resp[1,t_i,:]
			fields.hatR2[3,t_i,:] = 2*k3[2]*k3[3]*resp[1,t_i,:].*resp[2,t_i,:]
			
		end
		
	end
	
	return fields
	
end
