"""
	msfm()
	
	Calculate the shortest distance from a list of
	points to all other pixels in an image, using the
	Multistencil Fast Marching Method (MSFM). 
	
	This method gives more accurate
	distances by using second order derivatives and cross neighbours.
	
	T=msfm2d(F, SourcePoints, UseSecond, UseCross)
	
	# inputs,
	  F: The speed image. The speed function must always be larger
				than zero (min value 1e-8), otherwise some regions will
				never be reached because the time will go to infinity.
	  SourcePoints : A list of starting points [2 x N] (distance zero)
	  UseSecond : Boolean Set to true if not only first but also second
	               order derivatives are used (default)
	  UseCross : Boolean Set to true if also cross neighbours
	               are used (default)
	# outputs
	  T : Image with distance from SourcePoints to all pixels	
"""
function msfm(speedimage::AbstractArray{T,2}, SourcePointsIn::AbstractArray{T}, usesecond::Bool=true, usecross::Bool=true,Ed=false::Bool) where T

	distanceimage = zeros(T, size(speedimage))
	fill!(distanceimage, T(-Inf))
	# Augmented Fast Marching (For skeletonize)
	
	# Euclidian distance image
	if(Ed)
		Y = zeros(T,size(speedimage))
	end

	# Pixels which are processed and have a final distance are frozen
	Frozen = falses(size(speedimage));

	# Free memory to store neighbours of the (segmented) region
	neg_free = length(speedimage)
	neg_pos = 0
    neg_list_1 = zeros(T,neg_free)
    neg_list_x = fill(1, neg_free) 
    neg_list_y = fill(1, neg_free) 
	if Ed
        neg_list_4 = zeros(T,neg_free)
	end

	# (There are 3 pixel classes:
	#   - frozen (processed)
	#   - narrow band (boundary) (in list to check for the next pixel with smallest distance)
	#   - far (not yet used)

 	# Neighbours
	ne = [-1  0;
		   1  0;
		   0 -1;
		   0  1]

	# SourcePoints = Int.(floor.(SourcePointsIn))
	SourcePoints = SourcePointsIn

	# set all starting points to distance zero and frozen
	for z = 1:size(SourcePoints,2)
		# starting point
		x = Int(round(SourcePoints[1,z]))
		y = Int(round(SourcePoints[2,z]))
		# Set starting point to frozen and distance to zero
		Frozen[x,y] = true
		distanceimage[x,y] = 0
	end


	# Add all neighbours of the starting points to narrow list
	@inbounds for z = 1:size(SourcePoints,2)
		# starting point
		x = Int(round(SourcePoints[1,z]))
		y = Int(round(SourcePoints[2,z]))
		for k=1:4
			# Location of neighbour
			i = x + ne[k,1]
			j = y + ne[k,2]
			# Check if current neighbour is not yet frozen and inside the
			# picture
			if( (i>0) && (j>0) && (i<=size(speedimage,1)) && (j<=size(speedimage,2)) && (~Frozen[i,j]))
				Tt = 1/max(speedimage[i,j],eps(T))
				Ty = T(1)
				# Update distance in neighbour list or add to neighbour list
				if distanceimage[i,j] > 0
					if neg_list_1[distanceimage[i,j]] > Tt
						neg_list_1[distanceimage[i,j]] = Tt
					end
					if Ed
						neg_list_4[distanceimage[i,j]] = min(Ty, neg_list_4[distanceimage[i,j]]);
					end
				else
					neg_pos = neg_pos + 1
                    neg_list_1[neg_pos] = Tt
                    neg_list_x[neg_pos] = i
                    neg_list_y[neg_pos] = j
					if Ed
                        neg_list_4[neg_pos] = Ty
                    end
					distanceimage[i,j] = neg_pos
				end
			end
		end
	end

	# Reuse variables:
	Tpatch = zeros(T, 5, 5)
	Order = zeros(Int8, 1, 4)
	Tm = zeros(T,1,4)
	Tm2 = zeros(T,1,4)
	Coeff = zeros(T,3)
	TT = zeros(T,2)
	TT2 = zeros(T,2)
	
	# Loop through all pixels of the image
	for itt = 1:length(speedimage)
		# Get the pixel from narrow list (boundary list) with smallest
		# distance value and set it to current pixel location

		if neg_pos == 0
			break
		end
		(t,index) = findmin(view(neg_list_1,1:neg_pos))
		x = neg_list_x[index]
		y = neg_list_y[index]
		Frozen[x, y] = true
		distanceimage[x,y] = neg_list_1[index]

		if Ed
			Y[x,y] = neg_list_4[index]
		end

		# Remove min value by replacing it with the last value in the array
		if index < neg_pos
            neg_list_1[index] = neg_list_1[neg_pos]
            neg_list_x[index] = neg_list_x[neg_pos]
            neg_list_y[index] = neg_list_y[neg_pos]
            if Ed
                neg_list_4[index] = neg_list_4[neg_pos]
            end
			x2 = neg_list_x[index]
			y2 = neg_list_y[index]
			distanceimage[x2,y2] = index
		end
		neg_pos = neg_pos - 1
		# Loop through all 4 neighbours of current pixel
		@inbounds for k = 1:4
			# Location of neighbour
			i = x + ne[k,1]
			j = y + ne[k,2]

			# Check if current neighbour is not yet frozen and inside the
			# picture
			if (i>0) && (j>0) && (i<=size(speedimage,1)) && (j<= size(speedimage,2)) && (~Frozen[i,j])

				Tt = calculatedistance!(distanceimage, Tpatch, Order, Coeff, Tm, Tm2, TT, TT2,speedimage[i,j], size(speedimage), i, j, usesecond, usecross, Frozen)
				if Ed
					Ty = calculatedistance!(Y, Tpatch, Order, Coeff, Tm, Tm2, TT, TT2, speedimage[i,j], size(speedimage), i, j, usesecond, usecross, Frozen)
				end

				# Update distance in neighbour list or add to neighbour list
				if distanceimage[i,j] > 0
					if Tt !== nothing
						neg_list_1[round(Int,distanceimage[i,j])] = min(Tt, neg_list_1[round(Int,distanceimage[i,j])] )
					end
					if Ed && (Tt !== nothing)
						neg_list_4[round(Int,distanceimage[i,j])] = min(Ty,neg_list_4[round(Int,distanceimage[i,j])])
					end
				else
					neg_pos = neg_pos + 1
					neg_list_1[neg_pos] = Tt
					neg_list_x[neg_pos] = i
					neg_list_y[neg_pos] = j
					if Ed
						neg_list_4[neg_pos] = Ty
					end
					distanceimage[i,j] = neg_pos
				end
			end
		end
	end
	return distanceimage
end # End function msfm

function calculatedistance(TI::AbstractArray{T}, Fij, sizeF, i, j, usesecond, usecross, Frozen) where T
	Tpatch = zeros(T, 5, 5)
	Order = zeros(Int8, 1, 4)
	Tm = zeros(T,1,4)
	Tm2 = zeros(T,1,4)
	Coeff = zeros(T,3)
	TT = zeros(T,2)
	TT2 = zeros(T,2)
	calculatedistance!(TI, Tpatch, Order, Coeff, Tm, Tm2, TT, TT2, Fij, sizeF, i, j, usesecond, usecross, Frozen)
end


function calculatedistance!(TI::AbstractArray{T}, Tpatch, Order, Coeff, Tm, Tm2, TT, TT2, Fij, sizeF, i, j, usesecond, usecross, Frozen) where T
	
	fill!(Tpatch, T(Inf))
	@inbounds for nx = -2:2
		for ny = -2:2
			i_n = i + nx
			j_n = j + ny
			if i_n > 0 && j_n > 0 && i_n <= sizeF[1] && j_n <= sizeF[2] && Frozen[i_n, j_n]
				Tpatch[nx+3, ny+3] = TI[i_n, j_n]
			end
		end
	end

	# The values in order is 0 if no neighbours in that direction
	# 1 if 1e order derivatives is used and 2 if second order
	# derivatives are used
	fill!(Order, 0)

	# Make 1e order derivatives in x and y direction
	fill!(Tm, 0)
	Tm[1] = min(Tpatch[2,3], Tpatch[4,3] )
	if isfinite(Tm[1])
		Order[1]=1
	end

	Tm[2] = min(Tpatch[3,2], Tpatch[3,4])
	if isfinite(Tm[2])	
		Order[2]=1
	end

	# Make 1e order derivatives in cross directions
	if usecross
		Tm[3] = min(Tpatch[2,2], Tpatch[4,4])
		if isfinite(Tm[3])
			Order[3]=1
		end
		Tm[4] = min(Tpatch[2,4], Tpatch[4,2])
		if isfinite(Tm[4])
			Order[4]=1
		end
	end

	# Make 2e order derivatives
	if usesecond
		fill!(Tm2,0)
		# pixels with a pixeldistance 2 from the center must be
		# lower in value otherwise use other side or first order
		ch1 = (Tpatch[1,3] < Tpatch[2,3]) && isfinite(Tpatch[2,3])
		ch2 = (Tpatch[5,3] < Tpatch[4,3]) && isfinite(Tpatch[4,3])

		if ch1 && ch2
			Tm2[1] = min( (4*Tpatch[2,3] - Tpatch[1,3])/3 , (4*Tpatch[4,3] - Tpatch[5,3])/3 )
			Order[1] = 2
		elseif ch1
			Tm2[1] = (4*Tpatch[2,3] - Tpatch[1,3])/3
			Order[1] = 2
		elseif ch2
			Tm2[1] =(4*Tpatch[4,3] - Tpatch[5,3])/3
			Order[1] = 2
		end

		ch1 = (Tpatch[3,1] < Tpatch[3,2]) && isfinite(Tpatch[3,2])
		ch2 = (Tpatch[3,5] < Tpatch[3,4]) && isfinite(Tpatch[3,4])

		if ch1 && ch2
			Tm2[2] = min( (4*Tpatch[3,2] - Tpatch[3,1])/3 , (4*Tpatch[3,4] - Tpatch[3,5])/3)
			Order[2] = 2
		elseif ch1
			Tm2[2] = (4*Tpatch[3,2] - Tpatch[3,1])/3
			Order[2] = 2
		elseif ch2
			Tm2[2] = (4*Tpatch[3,4] - Tpatch[3,5])/3
			Order[2] = 2
		end

		if usecross
			ch1 = (Tpatch[1,1] < Tpatch[2,2]) && isfinite(Tpatch[2,2])
			ch2 = (Tpatch[5,5] < Tpatch[4,4]) && isfinite(Tpatch[4,4])
			if ch1 && ch2
				Tm2[3] = min( (4*Tpatch[2,2] - Tpatch[1,1])/3, (4*Tpatch[4,4] - Tpatch[5,5])/3)
				Order[3] = 2
			elseif ch1
				Tm2[3] = (4*Tpatch[2,2] - Tpatch[1,1])/3
				Order[3] = 2
			elseif ch2
				Tm2[3] = (4*Tpatch[4,4] - Tpatch[5,5])/3
				Order[3] = 2
			end

			ch1 = (Tpatch[1,5] < Tpatch[2,4]) && isfinite(Tpatch[2,4])
			ch2 = (Tpatch[5,1] < Tpatch[4,2]) && isfinite(Tpatch[4,2])
			if ch1 && ch2
				Tm2[4] = min( (4*Tpatch[2,4] - Tpatch[1,5])/3, (4*Tpatch[4,2] - Tpatch[5,1])/3)
				Order[4] = 2
			elseif ch1
				Tm2[4] = (4*Tpatch[2,4] - Tpatch[1,5])/3
				Order[4] = 2
			elseif ch2
				Tm2[4] = (4*Tpatch[4,2] - Tpatch[5,1])/3
				Order[4] = 2
			end
		end
	else
		fill!(Tm2,0)
	end

	# Calculate the distance using x and y direction
	Coeff .= (0, 0, -1/(max(Fij^2,eps(T))))
	for t = 1:2
		if Order[t] == 1
			Coeff .+= (1, -2*Tm[t], Tm[t]^2)
		elseif Order[t] == 2
			Coeff .+= (1, -2*Tm2[t], Tm2[t]^2) .* 2.25
		end
	end

	roots!(TT, Coeff)
	# Calculate the distance using the cross directions
	if usecross
		Coeff .+= (0, 0, -1/(maximum((Fij^2,eps(T)))))
		for t = 3:4
			if Order[t] == 1
				Coeff .+= 0.5 .* (1, -2*Tm[t], Tm[t]^2)
			elseif Order[t] == 2
				Coeff .+= 0.5 .* (1, -2*Tm2[t], Tm2[t]^2) .* 2.25
			end
		end
		roots!(TT2, Coeff)
		# Select minimum distance value of both stensils
		if ~isempty(TT2)
			TT = min(maximum(TT), maximum(TT2))
		end
	else
		TT = maximum(TT)
	end

	# Upwind condition check, current distance must be larger
	# then direct neighbours used in solution
	#DirectNeigbInSol=Tm(isfinite(Tm));
	#if(nnz(DirectNeigbInSol>=Tt)>0) # Will this ever happen?
	#    Tt=min(DirectNeigbInSol)+(1/(max(Fij,eps)));
	#end
end

function roots!(z, Coeff::AbstractArray{T}) where T
	a = Coeff[1]
	b = Coeff[2]
	c = Coeff[3]
	d = max((b*b)-4.0*a*c,zero(T))
	if a != 0
		z[1] = (-b - sqrt(d)) / (T(2.0)*a)
		z[2] = (-b + sqrt(d)) / (T(2.0)*a)
	else
		z[1] = (T(2.0)*c)/(-b - sqrt(d))
		z[2] = (T(2.0)*c)/(-b + sqrt(d))
	end
end

function roots_exp!(z, Coeff::AbstractArray{T}) where T
	a = Coeff[1]
	b = Coeff[2]
	c = Coeff[3]
	d = max((b*b)-4.0*a*c,zero(T))
	if a != 0
		z = max( (-b - sqrt(d)) / (T(2.0)*a), (-b + sqrt(d)) / (T(2.0)*a) )
	else
		z = max( (T(2.0)*c)/(-b - sqrt(d)), (T(2.0)*c)/(-b + sqrt(d)) )
	end
end

function roots(Coeff::AbstractArray{T}) where T
	z = zeros(T,2,1)
	roots!(z,Coeff)
	return z
end