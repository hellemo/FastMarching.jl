# RK4 is a function which performs one step of the Runge-Kutta 4 ray tracing
#
# EndPoint = RK4(StartPoint, GradientVolume, StepSize);
#
# inputs :
#      StartPoint: 2D or 3D location in vectorfield
#      GradientVolume: Vectorfield
#      Stepsize : The stepsize
#
# outputs :
#      EndPoint : The new location (zero if outside image)
#
# Function is written by D.Kroon University of Twente (July 2008)
# Port from C to Julia by Lars Hellemo (May 2014)

function checkBounds2d( point, Isize)
    if (point[1] < 1) || (point[2] < 1) || (point[1] > (Isize[1])) || (point[2] > (Isize[2]))
        return false
    end
    return true
end

function interpgrad2d!(Ireturn::AbstractArray{T},I, Isize, point) where T
    #  Linear interpolation variables
    perc = zeros(T,4,)
    index = fill(4,4,1)
    fTlocalx = floor(point[1])
    fTlocaly = floor(point[2])
    xBas0 = Int(fTlocalx)
    yBas0 = Int(fTlocaly)
    xBas1 = xBas0 + 1
    yBas1 = yBas0 + 1

    # Linear interpolation constants (percentages)
    xCom = point[1] - fTlocalx
    yCom = point[2] - fTlocaly
    xComi = (1 - xCom)
    yComi = (1 - yCom)
    perc[1] = xComi * yComi
    perc[2] = xComi * yCom
    perc[3] = xCom * yComi
    perc[4] = xCom * yCom

    # Stick to boundary
    if xBas0 < 1
        xBas0 = 1
        if xBas1 < 1
            xBas1 = 1
        end
    end

    if yBas0 < 1
        yBas0 = 1
        if yBas1 < 1
            yBas1 = 1
        end
    end
    if xBas1 > Isize[1]
        xBas1 = Isize[1]
        if xBas0 > Isize[1]
            xBas0 = Isize[1]
        end
    end

    if yBas1 > Isize[2]
        yBas1 = Isize[2]
        if yBas0 > Isize[2]
            yBas0 = Isize[2]
        end
    end

    # Get the neighbour intensities
    index[1] = LinearIndices(size(I))[xBas0,yBas0,1] #TODO Why is this dim 3,3,2?
    index[2] = LinearIndices(size(I))[xBas0,yBas1,1]
    index[3] = LinearIndices(size(I))[xBas1,yBas0,1]
    index[4] = LinearIndices(size(I))[xBas1,yBas1,1]

    f = Isize[1] * Isize[2]

    # the interpolated color
    Ireturn[1] = sum(I[index[i]] * perc[i] for i=1:4)
    Ireturn[2] = sum(I[index[i] + f] * perc[i] for i=1:4) 
    return Ireturn
end

"""
    Perform the RK4 raytracing step
"""
function rk4(startpoint::AbstractArray{T}, gradientvolume::AbstractArray{T}, stepsize::Number) where T
    gradientarray = gradientvolume # TODO Test if this works (without copy)
    gradientarraysize = collect(size(gradientarray))
    nextpoint = startpoint
    nextPoint =  rk4step2d(gradientarray, startpoint, ones(T,2,), stepsize)
    return nextPoint
end

"""
Deprecated name for rk4step2d
    """
function RK4STEP_2D(gradientArray::AbstractArray{T}, startPoint::AbstractArray{T}, nextPoint,stepSize) where T
    rk4step2d(gradientArray, startPoint, nextPoint, stepSize)
end
@deprecate RK4STEP_2D(gradientArray, startPoint, nextPoint, stepSize) rk4step2d(gradientarray, startpoint, nextpoint, stepsize)
"""
    Perform one step of the RK4 algorithm
"""
function rk4step2d(gradientarray::AbstractArray{T}, startpoint::AbstractArray{T}, nextpoint, stepsize) where T
    k = ones(T,4,2)
    temppoint = ones(T,2)
    gradientarraysize = collect(size(gradientarray))

    # Calculate k1:k4
    for ki = 1:4
        @views interpgrad2d!(k[ki,:], gradientarray, gradientarraysize, startpoint)
        @views tempnorm = max(norm(k[ki,:]),1e-6)
        for p = 1:2
            k[ki,p] *= stepsize/tempnorm
            temppoint[p] = startpoint[p] - k[ki,p] * T(0.5)
        end
        # Check if still inside the domain
        if !checkBounds2d(temppoint, gradientarraysize)
            return ones(T,2,)
        end
    end
    
    # Calculate final point
    for p = 1:2
        nextpoint[p] = startpoint[p] - (k[1,p] + k[2,p] * T(2.0) + k[3,p] * T(2.0) + k[4,p])/T(6.0)
    end

    # Set step to step size
    # D1 = nextPoint[1] - startPoint[1]
    # D2 = nextPoint[2] - startPoint[2]
    # dl = stepSize/(sqrt(D1 * D1 + D2 * D2)+1e-15)
    # D1 *= dl
    # D2 *= dl
    # nextPoint[1] = startPoint[1] + D1
    # nextPoint[2] = startPoint[2] + D2

    # Check the if are still inside the domain
    if !checkBounds2d(nextpoint, gradientarraysize)
        return ones(T,2,)
    end
    return nextpoint
end
