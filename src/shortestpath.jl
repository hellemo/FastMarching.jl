"""
    shortestpath(distancemap, startpoint, sourcepoint, stepsize)

    This function traces the shortest path from `startpoint` to
    `sourcepoint` using Runge Kutta 4 in a 2D `distancemap`.

    Optionally specify the line tracing `stepsize`.

    Returns a 2D-array with the shortest path.
"""
function shortestpath(distancemap::AbstractArray{T}, startpoint::AbstractArray{T}, sourcepoint=ones(T,2,1), stepsize::T=T(0.5)) where T
    gradientvolume = zeros(T, size(distancemap,1), size(distancemap,2),2)
    shortestpath!(gradientvolume, distancemap, startpoint, sourcepoint, stepsize)
end

function shortestpath!(gradientvolume, distancemap::AbstractArray{T}, startpoint::AbstractArray{T}, sourcepoint=ones(T,2,1), stepsize::T=0.5) where T
    
    distancetoend = T(Inf)    
    Fx = view(gradientvolume,:,:,1)
    Fy = view(gradientvolume,:,:,2)
    pointmin!(distancemap, Fx, Fy)
    Fx .= -Fx
    Fy .= -Fy
    
    i = 0
    # Reserve a block of memory for the shortest line array
    ifree = size(distancemap, 1) * 2
    shortestline = zeros(T, ifree, ndims(distancemap))

    # Iteratively trace the shortest line
    movement = T(1)
    while true
        endpoint = rk4(startpoint, gradientvolume, stepsize)

        # Calculate the distance to the end point
        if !isempty(sourcepoint) && !isnan(sourcepoint[1])
            (distancetoend,ind) = findmin(sqrt.(sum((sourcepoint.-repeat(endpoint,1,size(sourcepoint,2))).^2)))
        else
            distancetoend = Inf
        end

        # Calculate the movement between current point and point 10 iterations back
        if i > 10
            movement = sqrt(sum((endpoint - shortestline[i-10,:]).^2))
        end

        # Stop if out of boundary, distance to end smaller then a pixel or
        # if we have not moved for 10 iterations
        if (endpoint[1] == 0) || (movement < stepsize)
            break
        end

        # Count the number of iterations
        i = i + 1

        # Add a new block of memory if nearly full
        if i > (ifree - 2)
            ifree = ifree * 2
            shortestline = vcat(shortestline,zeros(T,ifree, ndims(distancemap))) 
        end

        # Add current point to the shortest line array
        shortestline[i,:] = endpoint

        if distancetoend < stepsize
            i = i + 1
            # Add (Last) Source point to the shortest line array
            shortestline[i,:] = sourcepoint
            break
        end

        # Current point is next Starting Point
        startpoint = endpoint
    end

    if (distancetoend  > 1) && (~isempty(sourcepoint))
        @warn "The shortest path trace did not finish at the source point" distancetoend
    end

    # Remove unused memory from array
    shortestline = shortestline[1:i,:]
    return shortestline
end
