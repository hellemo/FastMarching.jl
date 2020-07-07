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
function msfm(speedimage::AbstractArray{T,2}, SourcePoints::AbstractArray{T}, usesecond::Bool=true, usecross::Bool=true,Ed=false::Bool) where T

	distanceimage = fill(T(-1.0),size(speedimage))
	# Augmented Fast Marching (For skeletonize)
	#Ed=false # Was: nargout>1;

	# Euclidian distance image
	if(Ed)
		Y = zeros(T,size(speedimage))
	end

	# Pixels which are processed and have a final distance are frozen
	Frozen   = falses(size(speedimage));

	# Free memory to store neighbours of the (segmented) region
	neg_free = 100000
	neg_pos=0
	if Ed
		neg_list = zeros(4,neg_free)
	else
		neg_list = zeros(3,neg_free)
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

	SourcePoints=round.(Int,floor.(SourcePoints))


	# set all starting points to distance zero and frozen
	for z=1:size(SourcePoints,2)
		# starting point
		x= SourcePoints[1,z]
		y= SourcePoints[2,z]
		# Set starting point to frozen and distance to zero
		Frozen[x,y]=true
		distanceimage[x,y]=zero(T)
	end


	# Add all neighbours of the starting points to narrow list
	for z = 1:size(SourcePoints,2)
		# starting point
		x = SourcePoints[1,z]
		y = SourcePoints[2,z]
		for k=1:4
			# Location of neighbour
			i = x + ne[k,1]
			j = y + ne[k,2]
			# Check if current neighbour is not yet frozen and inside the
			# picture
			if( (i>0) && (j>0) && (i<=size(speedimage,1)) && (j<=size(speedimage,2)) && (~Frozen[i,j]))
				Tt = 1/max(speedimage[i,j],eps(T))
				Ty = T(1)
				# Update distance in neigbour list or add to neigbour list
				if distanceimage[i,j] > 0
					if neg_list[1,distanceimage[i,j]] > Tt
						neg_list[1,distanceimage[i,j]] = Tt
					end
					if Ed
						neg_list[4,distanceimage[i,j]] = min(Ty, neg_list[4,distanceimage[i,j]]);
					end
				else
					neg_pos=neg_pos+1
					# If running out of memory at a new block
					if neg_pos>neg_free
						neg_free = neg_free + 100_000
						if Ed
							neg_list = hcat(neg_list, zeros(4,neg_free))
						else
							neg_list = hcat(neg_list, zeros(3,neg_free))
						end
					end
					if Ed
						@views neg_list[:,neg_pos] = [Tt; i; j; Ty]
					else
						@views neg_list[:,neg_pos] = [Tt; i; j]
					end
					distanceimage[i,j] = neg_pos
				end
			end
		end
	end


	# Reuse variables for all iterations:
	Tpatch = @MArray zeros(T, 5, 5)
	Order = @MArray zeros(Int8, 1, 4)
	Tm = @MArray zeros(T,1,4)
	Tm2 = @MArray zeros(T,1,4)
	Coeff = @MArray zeros(T,3)

	# Loop through all pixels of the image
	for itt = 1:length(speedimage)
		# Get the pixel from narrow list (boundary list) with smallest
		# distance value and set it to current pixel location

		if neg_pos==0
			break
		end
		(t,index) = findmin(neg_list[1,1:neg_pos])
		x = round(Int,neg_list[2,index])
		y = round(Int,neg_list[3,index])
		Frozen[x, y] = true
		distanceimage[x,y] = neg_list[1,index]

		if Ed
			Y[x,y] = neg_list[4,index]
		end

		# Remove min value by replacing it with the last value in the array
		if index < neg_pos
			@views neg_list[:,index] = neg_list[:,neg_pos]
			x2 = round(Int,neg_list[2,index])
			y2 = round(Int,neg_list[3,index])
			distanceimage[x2,y2] = index
		end
		neg_pos = neg_pos-1

		# Loop through all 4 neighbours of current pixel
		for k = 1:4
			# Location of neighbour
			i = x + ne[k,1]
			j = y + ne[k,2]

			# Check if current neighbour is not yet frozen and inside the
			# picture
			if (i>0)&&(j>0)&&(i<=size(speedimage,1))&&(j<=size(speedimage,2))&&(~Frozen[i,j])

				Tt = calculatedistance!(distanceimage, Tpatch, Order, Coeff, Tm, Tm2, speedimage[i,j], size(speedimage), i, j, usesecond, usecross, Frozen)
				if Ed
					Ty = calculatedistance!(Y, Tpatch, Order, Coeff, Tm, Tm2, speedimage[i,j], size(speedimage), i, j, usesecond, usecross, Frozen)
				end

				# Update distance in neigbour list or add to neigbour list
				if distanceimage[i,j]>0
					neg_list[1,round(Int,distanceimage[i,j])] = minimum((Tt,neg_list[1,round(Int,distanceimage[i,j])]))
					if Ed
						neg_list[4,round(Int,distanceimage[i,j])] = minimum((Ty,neg_list[4,round(Int,distanceimage[i,j])]))
					end
				else
					neg_pos = neg_pos + 1
					# If running out of memory at a new block
					if neg_pos>neg_free
						neg_free = neg_free +100_000
						#Was: neg_list(1,neg_free)=0;
						if Ed
							neg_list = hcat(neg_list,zeros(4,neg_free))
						else
							neg_list = hcat(neg_list,zeros(3,neg_free))
						end
					end
					if Ed
						@views neg_list[:,neg_pos] = [Tt; i; j; Ty]
					else
						@views neg_list[:,neg_pos]=[Tt; i; j]
					end
					distanceimage[i,j] = neg_pos
				end
			end
		end
	end
	return distanceimage
end # End function msfm

function calculatedistance(TI::AbstractArray{T}, Fij, sizeF, i, j, usesecond, usecross, Frozen) where T
	Tpatch = @MArray zeros(T, 5, 5)
	Order = @MArray zeros(Int8, 1, 4)
	Tm = @MArray zeros(T,1,4)
	Tm2 = @MArray zeros(T,1,4)
	Coeff = @MArray zeros(T,3)
	# Tt = zero(T)
	# Tt2 = zero(T)
	calculatedistance!(TI, Tpatch, Order, Coeff, Tm, Tm2, Fij, sizeF, i, j, usesecond, usecross, Frozen)
end


function calculatedistance!(TI::AbstractArray{T}, Tpatch, Order, Coeff, Tm, Tm2, Fij, sizeF, i, j, usesecond, usecross, Frozen) where T
	
	fill!(Tpatch,T(Inf))
	for nx = -2:2
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
	fill!(Order,0)

	# Make 1e order derivatives in x and y direction
	fill!(Tm,0)
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
			Tm2[1] = minimum(( (4*Tpatch[2,3] - Tpatch[1,3])/3 , (4*Tpatch[4,3] - Tpatch[5,3])/3 ))
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
				Tm2[3] = min( (4*Tpatch[2,2] - Tpatch[1,1])/3 , (4*Tpatch[4,4] - Tpatch[5,5])/3)
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
				Tm2[4] = minimum(( (4*Tpatch[2,4] - Tpatch[1,5])/3 , (4*Tpatch[4,2] - Tpatch[5,1])/3))
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
		Tm2=zeros(T,1,4)
	end

	# Calculate the distance using x and y direction
	# Coeff = [0 0 -1/(max(Fij^2,eps()))]
	Coeff .= T.([0, 0, -1/(max(Fij^2,eps(T)))])
	for t = 1:2
		if Order[t] == 1
			Coeff .+= [1, -2*Tm[t], Tm[t]^2]
		elseif Order[t] == 2
			Coeff .+= [1, -2*Tm2[t], Tm2[t]^2] .* 2.25
		end
	end

	Tt=roots(Coeff)
	Tt=maximum(Tt)
	# Calculate the distance using the cross directions
	if usecross
		Coeff .+= T.([0, 0, -1/(maximum([Fij^2,eps(T)]))])
		for t=3:4
			if Order[t] == 1
				Coeff .+= 0.5 .*[1, -2*Tm[t], Tm[t]^2]
			elseif Order[t] == 2
				Coeff .+= 0.5 .*[1, -2*Tm2[t], Tm2[t]^2] .* 2.25
			end
		end
		Tt2=roots(Coeff)
		Tt2=maximum(Tt2)
		# Select minimum distance value of both stensils
		if ~isempty(Tt2)
			Tt=minimum((Tt,Tt2))
		end
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
	z = @MArray zeros(T,2,1)
	roots!(z,Coeff)
	return z
end