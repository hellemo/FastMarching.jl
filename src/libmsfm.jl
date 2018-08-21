# TODO: Check usefulness/relevant size of neg_free
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
	
	# Literature
	  M. Sabry Hassouna et Al. Multistencils Fast Marching
	  Methods: A Highly Accurate Solution to the Eikonal Equation on
	  Cartesian Domains
	
	# Example
	  SourcePoint = [51; 51];
	  SpeedImage = ones([101 101]);
	  [X Y] = ndgrid(1:101, 1:101);
	  T1 = sqrt((X-SourcePoint(1)).^2 + (Y-SourcePoint(2)).^2);
	
	  # Run fast marching 1th order, 1th order multi stencil
	  # and 2th orde and 2th orde multi stencil
	
	  tic; T1_FMM1 = msfm2d(SpeedImage, SourcePoint, false, false); toc;
	  tic; T1_MSFM1 = msfm2d(SpeedImage, SourcePoint, false, true); toc;
	  tic; T1_FMM2 = msfm2d(SpeedImage, SourcePoint, true, false); toc;
	  tic; T1_MSFM2 = msfm2d(SpeedImage, SourcePoint, true, true); toc;
	
	  # Show results
	  fprintf('\nResults with T1 (Matlab)\n');
	  fprintf('Method   L1        L2        Linf\n');
	  Results = cellfun(@(x)([mean(abs(T1(:)-x(:))) mean((T1(:)-x(:)).^2) max(abs(T1(:)-x(:)))]), {T1_FMM1(:) T1_MSFM1(:) T1_FMM2(:) T1_MSFM2(:)}, 'UniformOutput',false);
	  fprintf('FMM1:   %9.5f %9.5f %9.5f\n', Results{1}(1), Results{1}(2), Results{1}(3));
	  fprintf('MSFM1:  %9.5f %9.5f %9.5f\n', Results{2}(1), Results{2}(2), Results{2}(3));
	  fprintf('FMM2:   %9.5f %9.5f %9.5f\n', Results{3}(1), Results{3}(2), Results{3}(3));
	  fprintf('MSFM2:  %9.5f %9.5f %9.5f\n', Results{4}(1), Results{4}(2), Results{4}(3));
	
	# Example multiple starting points,
	  SourcePoint=rand(2,100)*255+1;
	  SpeedImage = ones([256 256]);
	  tic; T1_MSFM2 = msfm2d(SpeedImage, SourcePoint, true, true); toc;
	  figure, imshow(T1_MSFM2,[]); colormap(hot(256));
	
	  # Credit
	  Function is written by D.Kroon University of Twente (June 2009)
	  Port from Matlab to Julia by Lars Hellemo (May 2014)
"""
function msfm(F::Array{Float64,2}, SourcePoints, usesecond::Bool, usecross::Bool,Ed=false::Bool)
	
	# Distance image, also used to store the index of narrowband pixels
	# during marching process
	T = fill(-1.0,size(F))
	# Augmented Fast Marching (For skeletonize)
	#Ed=false # Was: nargout>1;

	# Euclidian distance image
	if(Ed)
		Y = zeros(size(F))
	end

	# Pixels which are processed and have a final distance are frozen
	Frozen   = falses(size(F));

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
	ne =[-1 0;
	1 0;
	0 -1;
	0 1];

	SourcePoints=round.(Int,floor.(SourcePoints))


	# set all starting points to distance zero and frozen
	for z=1:size(SourcePoints,2)
		# starting point
		x= SourcePoints[1,z]
		y= SourcePoints[2,z]
		# Set starting point to frozen and distance to zero
		Frozen[x,y]=true
		T[x,y]=0
	end




	# Add all neighbours of the starting points to narrow list
	for z=1:size(SourcePoints,2)
		# starting point
		x=SourcePoints[1,z]
		y=SourcePoints[2,z]
		for k=1:4
			# Location of neighbour
			i=x+ne[k,1]
			j=y+ne[k,2]
			# Check if current neighbour is not yet frozen and inside the
			# picture
			if((i>0)&&(j>0)&&(i<=size(F,1))&&(j<=size(F,2))&&(~Frozen[i,j]))
				Tt=1/max(F[i,j],eps())
				Ty=1;
				# Update distance in neigbour list or add to neigbour list
				if T[i,j]>0
					if neg_list[1,T[i,j]]>Tt
						neg_list[1,T[i,j]]=Tt
					end
					if Ed
						neg_list[4,T[i,j]]=min(Ty,neg_list[4,T[i,j]]);
					end
				else
					neg_pos=neg_pos+1
					# If running out of memory at a new block
					if neg_pos>neg_free
						neg_free = neg_free + 100000
						if Ed
							neg_list = hcat(neg_list,zeros(4,neg_free))# Was: neg_list(1,neg_free)=0;
						else
							neg_list = hcat(neg_list,zeros(3,neg_free))# Was: neg_list(1,neg_free)=0;
						end
					end
					if Ed
						neg_list[:,neg_pos]=[Tt;i;j;Ty]
					else
						neg_list[:,neg_pos]=[Tt;i;j]
					end
					T[i,j]=neg_pos
				end
			end
		end
	end


	# Loop through all pixels of the image
	for itt=1:length(F)
		# Get the pixel from narrow list (boundary list) with smallest
		# distance value and set it to current pixel location

		if neg_pos==0
			break
		end
		(t,index)=findmin(neg_list[1,1:neg_pos])
		x=round(Int,neg_list[2,index])
		y=round(Int,neg_list[3,index])
		Frozen[x,y]=true
		# T[x,y]=round(Int,neg_list[1,index])
		T[x,y]=neg_list[1,index]

		if Ed
			Y[x,y]=neg_list[4,index]
		end

		# Remove min value by replacing it with the last value in the array
		if index<neg_pos
			neg_list[:,index]=neg_list[:,neg_pos]
			x2=round(Int,neg_list[2,index])
			y2=round(Int,neg_list[3,index])
			T[x2,y2]=index
		end
		neg_pos =neg_pos-1

		# Loop through all 4 neighbours of current pixel
		for k=1:4
			# Location of neighbour
			i=x+ne[k,1]
			j=y+ne[k,2]

			# Check if current neighbour is not yet frozen and inside the
			# picture
			if (i>0)&&(j>0)&&(i<=size(F,1))&&(j<=size(F,2))&&(~Frozen[i,j])

				Tt=CalculateDistance(T,F[i,j],size(F),i,j,usesecond,usecross,Frozen)
				if Ed
					Ty=CalculateDistance(Y,1,size(F),i,j,usesecond,usecross,Frozen)
				end

				# Update distance in neigbour list or add to neigbour list
				if T[i,j]>0
					neg_list[1,round(Int,T[i,j])]=minimum([Tt,neg_list[1,round(Int,T[i,j])]])
					if Ed
						neg_list[4,round(Int,T[i,j])]=minimum([Ty,neg_list[4,round(Int,T[i,j])]])
					end
				else
					neg_pos=neg_pos+1
					# If running out of memory at a new block
					if neg_pos>neg_free
						neg_free = neg_free +100_000
						#Was: neg_list(1,neg_free)=0;
						if Ed
							neg_list = hcat(neg_list,zeros(4,neg_free))# Was: neg_list(1,neg_free)=0;
						else
							neg_list = hcat(neg_list,zeros(3,neg_free))# Was: neg_list(1,neg_free)=0;
						end
					end
					if Ed
						neg_list[:,neg_pos]=[Tt;i;j;Ty]
					else
						neg_list[:,neg_pos]=[Tt;i;j]
					end
					T[i,j]=neg_pos
				end
			end
		end
	end

	return T

end # End function msfm







function CalculateDistance(T,Fij,sizeF,i,j,usesecond,usecross,Frozen)
	#  Tpatch=infs(5,5)
	Tpatch = fill(Inf,(5,5))
	for nx = -2:2
		for ny = -2:2
			i_n = i + nx
			j_n = j + ny
			if i_n>0 && j_n>0 && i_n<=sizeF[1] && j_n<=sizeF[2] && Frozen[i_n,j_n]
				Tpatch[nx+3,ny+3] = T[i_n,j_n]
			end
		end
	end


	# The values in order is 0 if no neighbours in that direction
	# 1 if 1e order derivatives is used and 2 if second order
	# derivatives are used
	Order=zeros(Int,1,4);


	# Make 1e order derivatives in x and y direction
	Tm=zeros(1,4)
	Tm[1] = min( Tpatch[2,3] , Tpatch[4,3] )
	if isfinite(Tm[1])
		Order[1]=1
	end

	Tm[2] = min( Tpatch[3,2] , Tpatch[3,4])
	if isfinite(Tm[2])
		Order[2]=1
	end

	# Make 1e order derivatives in cross directions
	if usecross
		Tm[3] = minimum( [Tpatch[2,2] , Tpatch[4,4]])
		if isfinite(Tm[3])
			Order[3]=1
		end
		Tm[4] = min( Tpatch[2,4] , Tpatch[4,2])
		if isfinite(Tm[4])
			Order[4]=1
		end
	end



	# Make 2e order derivatives
	if usesecond
		Tm2=zeros(1,4)
		# pixels with a pixeldistance 2 from the center must be
		# lower in value otherwise use other side or first order
		ch1 = (Tpatch[1,3]<Tpatch[2,3])&&isfinite(Tpatch[2,3])
		ch2 = (Tpatch[5,3]<Tpatch[4,3])&&isfinite(Tpatch[4,3])

		if ch1&&ch2
			Tm2[1] =minimum([ (4*Tpatch[2,3]-Tpatch[1,3])/3 , (4*Tpatch[4,3]-Tpatch[5,3])/3 ])
			Order[1]=2
		elseif ch1
			Tm2[1] =(4*Tpatch[2,3]-Tpatch[1,3])/3
			Order[1]=2
		elseif ch2
			Tm2[1] =(4*Tpatch[4,3]-Tpatch[5,3])/3
			Order[1] = 2
		end

		ch1 = (Tpatch[3,1]<Tpatch[3,2])&&isfinite(Tpatch[3,2])
		ch2 = (Tpatch[3,5]<Tpatch[3,4])&&isfinite(Tpatch[3,4])

		if ch1&&ch2
			Tm2[2] = min( (4*Tpatch[3,2]-Tpatch[3,1])/3 , (4*Tpatch[3,4]-Tpatch[3,5])/3)
			Order[2]=2
		elseif ch1
			Tm2[2]=(4*Tpatch[3,2]-Tpatch[3,1])/3
			Order[2]=2
		elseif ch2
			Tm2[2]=(4*Tpatch[3,4]-Tpatch[3,5])/3
			Order[2]=2
		end

		if usecross
			ch1 = (Tpatch[1,1]<Tpatch[2,2])&&isfinite(Tpatch[2,2])
			ch2 = (Tpatch[5,5]<Tpatch[4,4])&&isfinite(Tpatch[4,4])
			if ch1&&ch2
				Tm2[3] =min( (4*Tpatch[2,2]-Tpatch[1,1])/3 , (4*Tpatch[4,4]-Tpatch[5,5])/3)
				Order[3]=2
			elseif ch1
				Tm2[3]=(4*Tpatch[2,2]-Tpatch[1,1])/3
				Order[3]=2
			elseif ch2
				Tm2[3]=(4*Tpatch[4,4]-Tpatch[5,5])/3
				Order[3]=2
			end

			ch1 = (Tpatch[1,5]<Tpatch[2,4])&&isfinite(Tpatch[2,4])
			ch2 = (Tpatch[5,1]<Tpatch[4,2])&&isfinite(Tpatch[4,2])
			if ch1&&ch2
				Tm2[4] = minimum([ (4*Tpatch[2,4]-Tpatch[1,5])/3 , (4*Tpatch[4,2]-Tpatch[5,1])/3])
				Order[4]=2
			elseif ch1
				Tm2[4]=(4*Tpatch[2,4]-Tpatch[1,5])/3
				Order[4]=2
			elseif ch2
				Tm2[4] = (4*Tpatch[4,2]-Tpatch[5,1])/3
				Order[4]=2
			end
		end
	else
		Tm2=zeros(1,4)
	end


	# Calculate the distance using x and y direction
	Coeff = [0 0 -1/(max(Fij^2,eps()))]
	for t=1:2
		if Order[t] == 1
			Coeff=Coeff+[1 -2*Tm[t] Tm[t]^2]
		elseif Order[t] ==2
			Coeff=Coeff+[1 -2*Tm2[t] Tm2[t]^2]*(2.2500)
		end
	end


	Tt=roots(Coeff)
	Tt=maximum(Tt)
	# Calculate the distance using the cross directions
	if usecross
		Coeff = Coeff + [0 0 -1/(maximum([Fij^2,eps()]))]
		for t=3:4
			if Order[t] == 1
				Coeff=Coeff+0.5*[1 -2*Tm[t] Tm[t]^2]
			elseif Order[t] == 2
				Coeff=Coeff+0.5*[1 -2*Tm2[t] Tm2[t]^2]*(2.2500)
			end
		end
		Tt2=roots(Coeff)
		Tt2=maximum(Tt2)
		# Select minimum distance value of both stensils
		if ~isempty(Tt2)
			Tt=minimum([Tt,Tt2])
		end
	end

	# Upwind condition check, current distance must be larger
	# then direct neighbours used in solution
	#DirectNeigbInSol=Tm(isfinite(Tm));
	#if(nnz(DirectNeigbInSol>=Tt)>0) # Will this ever happen?
	#    Tt=min(DirectNeigbInSol)+(1/(max(Fij,eps)));
	#end


end

function roots(Coeff)
	a=Coeff[1]
	b=Coeff[2]
	c=Coeff[3]
	d=maximum([(b*b)-4.0*a*c,0])
	z = zeros(2,1)
	if a != 0
		z[1]= (-b - sqrt(d)) / (2.0*a)
		z[2]= (-b + sqrt(d)) / (2.0*a)
	else
		z[1]= (2.0*c)/(-b - sqrt(d))
		z[2]= (2.0*c)/(-b + sqrt(d))
	end
end


function run_test()
	#	 include("ndgrid.jl")

	# Example/test
	#SourcePoint = [51; 51]
	#SpeedImage = ones(101, 101)
	#XY = ndgrid(1:101, 1:101)
	#T1 = sqrt((XY[1]-SourcePoint[1]).^2 + (XY[2]-SourcePoint[2]).^2)
	#tic()
	#T1_MSFM2 = msfm2d(SpeedImage, SourcePoint, true, true)
	#toc()

	# Example multiple starting points,
	SourcePoint=rand(2,100)*255+1
	SpeedImage = ones(256, 256)
	tic()
	T1_MSFM2 = msfm(SpeedImage, SourcePoint, true, true)
	toc()
	writecsv("tmp.csv",T1_MSFM2)
end

#    run_test()
