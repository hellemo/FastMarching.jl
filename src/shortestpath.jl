
function shortestpath(DistanceMap::Array{Float64,2},StartPoint::Array{Float64,1},SourcePoint=ones(2,)::Array{Float64,1},Stepsize=0.5::Float64,Method="rk4")
    # This function SHORTESTPATH traces the shortest path from start point to
    # source point using Runge Kutta 4 in a 2D or 3D distance map.
    #
    # ShortestLine=shortestpath(DistanceMap,StartPoint,SourcePoint,Stepsize,Method)
    #
    # inputs,
    #   DistanceMap : A 2D or 3D distance map (from the functions msfm2d or msfm3d)
    #   StartPoint : Start point of the shortest path
    #   SourcePoint : (Optional), End point of the shortest path
    #   Stepsize: (Optional), Line trace step size
    #   Method: (Optional), 'rk4' (default), 'euler' ,'simple'
    # output,
    #   ShortestLine: M x 2 or M x 3 array with the Shortest Path
    #
    # Example,
    #   # Load a maze image
    #   I1=im2double(imread('images/maze.gif'));
    #   # Convert the image to a speed map
    #   SpeedImage=I1*1000+0.001;
    #   # Set the source to end of the maze
    #   SourcePoint=[800;803];
    #   # Calculate the distance map (distance to source)
    #   DistanceMap= msfm(SpeedImage, SourcePoint);
    #   # Show the distance map
    #   figure, imshow(DistanceMap,[0 3400])
    #   # Trace shortestline from StartPoint to SourcePoint
    #   StartPoint=[9;14];
    #   ShortestLine=shortestpath(DistanceMap,StartPoint,SourcePoint);
    #   # Plot the shortest route?
    DistancetoEnd = 1e300#Inf
    GradientVolume = zeros(size(DistanceMap,1),size(DistanceMap,2),2)
    #if ndims(DistanceMap) == 2 # Select 2D or 3D
        (Fx,Fy) = pointmin(DistanceMap)
        GradientVolume[:,:,1]=-Fx
        GradientVolume[:,:,2]=-Fy
    #else
    #    (Fy,Fx,Fz) = pointmin(DistanceMap)
    #    GradientVolume[:,:,:,1]=-Fx
    #    GradientVolume[:,:,:,2]=-Fy
    #    GradientVolume[:,:,:,3]=-Fz
    #end

    i=0
    # Reserve a block of memory for the shortest line array
    ifree=10000
    ShortestLine=zeros(ifree,ndims(DistanceMap))

    # Iteratively trace the shortest line
    while true
        # TODO: Implement other methods than Runge Kutta
        # # Calculate the next point using runge kutta
        # switch(lower(Method))
        #     case 'rk4'
        #         EndPoint=rk4(StartPoint, GradientVolume, Stepsize);
        #     case 'euler'
        #         EndPoint=e1(StartPoint, GradientVolume, Stepsize);
        #     case 'simple'
        #         EndPoint=s1(StartPoint,DistanceMap);
        #     otherwise
        #         error("shortestpath:input","unknown method");
        # end

#        EndPoint = s1(StartPoint,DistanceMap)
        EndPoint = rk4(StartPoint, GradientVolume, Stepsize)
#        EndPoint = e1(StartPoint, GradientVolume, Stepsize)

        # Calculate the distance to the end point
        if !isempty(SourcePoint) && !isnan(SourcePoint[1])
            (DistancetoEnd,ind)=findmin(sqrt.(sum((SourcePoint.-repeat(EndPoint,1,size(SourcePoint,2))).^2)))
        else
#            DistancetoEnd=Inf
            DistancetoEnd = Inf
        end

        # Calculate the movement between current point and point 10 iterations back
        if i>10
            Movement=sqrt(sum((EndPoint[:]-ShortestLine[i-10,:]).^2))
        else
            Movement=Stepsize+1
        end

        # Stop if out of boundary, distance to end smaller then a pixel or
        # if we have not moved for 10 iterations
        if (EndPoint[1]==0)||(Movement<Stepsize)
            break
        end

        # Count the number of iterations
        i=i+1

        # Add a new block of memory if full
        if i>ifree
            ifree=ifree+10000
            ShortestLine = vcat(ShortestLine,zeros(10000,ndims(DistanceMap))) # Was: ShortestLine(ifree,:)=0
        end

        # Add current point to the shortest line array
        ShortestLine[i,:]=EndPoint

        if DistancetoEnd<Stepsize
            i=i+1
            if i>ifree
                ifree=ifree+10000
                ShortestLine = vcat(ShortestLine,zeros(10000,ndims(DistanceMap))) # Was: ShortestLine(ifree,:)=0
            end
            # Add (Last) Source point to the shortest line array
            ShortestLine[i,:]=SourcePoint#[:,ind]
            break
        end

        # Current point is next Starting Point
        StartPoint=EndPoint
    end


    if (DistancetoEnd>1)&&(~isempty(SourcePoint))
        print("The shortest path trace did not finish at the source point");
    end

    # Remove unused memory from array
    ShortestLine=ShortestLine[1:i,:]
    return ShortestLine
end
