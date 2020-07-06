"""
    shortestpath(distancemap, startpoint, sourcepoint, stepsize)

    This function traces the shortest path from `startpoint` to
    `sourcepoint` using Runge Kutta 4 in a 2D `distancemap`.

    Optionally specify the line tracing `stepsize`.

    Returns a 2D-array with the shortest path.
"""
function shortestpath(distancemap::AbstractArray{T},startpoint::AbstractArray{T},sourcepoint=ones(T,2,1),stepsize::T=T(0.5)) where T
    
    distancetoend = T(Inf)
    gradientvolume = zeros(T,size(distancemap,1),size(distancemap,2),2)
    (Fx,Fy) = pointmin(distancemap)
    @views gradientvolume[:,:,1]=-Fx
    @views gradientvolume[:,:,2]=-Fy
    
    i=0
    # Reserve a block of memory for the shortest line array
    ifree=size(distancemap,1) * 2
    ShortestLine=zeros(T,ifree,ndims(distancemap))

    # Iteratively trace the shortest line
    while true
        EndPoint = rk4(startpoint, gradientvolume, stepsize)

        # Calculate the distance to the end point
        if !isempty(sourcepoint) && !isnan(sourcepoint[1])
            (distancetoend,ind)=findmin(sqrt.(sum((sourcepoint.-repeat(EndPoint,1,size(sourcepoint,2))).^2)))
        else
            distancetoend = Inf
        end

        # Calculate the movement between current point and point 10 iterations back
        if i>10
            Movement=sqrt(sum((EndPoint[:]-ShortestLine[i-10,:]).^2))
        else
            Movement=stepsize+1
        end

        # Stop if out of boundary, distance to end smaller then a pixel or
        # if we have not moved for 10 iterations
        if (EndPoint[1]==0)||(Movement<stepsize)
            break
        end

        # Count the number of iterations
        i=i+1

        # Add a new block of memory if nearly full
        if i > (ifree - 2)
            ifree=ifree * 2
            ShortestLine = vcat(ShortestLine,zeros(T,ifree, ndims(distancemap))) 
        end

        # Add current point to the shortest line array
        ShortestLine[i,:]=EndPoint

        if distancetoend<stepsize
            i=i+1
            # Add (Last) Source point to the shortest line array
            ShortestLine[i,:] = sourcepoint
            break
        end

        # Current point is next Starting Point
        startpoint=EndPoint
    end

    if (distancetoend>1)&&(~isempty(sourcepoint))
        @warn "The shortest path trace did not finish at the source point"
    end

    # Remove unused memory from array
    ShortestLine=ShortestLine[1:i,:]
    return ShortestLine
end
