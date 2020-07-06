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


#function mindex2(x, y, sizx)
#    return y*sizx+x
#end

function mindex3(x::Int, y::Int, z::Int, sizx::Int, sizy::Int)
    return z*sizy*sizx+y*sizx+x
end

function checkBounds2d( point, Isize)
    if (point[1]<1)||(point[2]<1)||(point[1]>(Isize[1]))||(point[2]>(Isize[2]))
        return false
    end
    return true
end

# TODO check dimensions of point, seems to be bogus
function checkBounds3d( point::Array{Number,3}, Isize::Array{Int,1})
    if (point[1,1,1]<1)||(point[2,1,1]<1)||(point[3,1,1]<1)||(point[1,1,1]>(Isize[1]))||(point[2,1,1]>(Isize[2]))||(point[3,1,1]>(Isize[3]))
         return false
    end
    return true
end

function interpgrad2d(Ireturn::AbstractArray{T},I, Isize, point) where T
    #  Linear interpolation variables

    perc = zeros(T,4,)
    index = fill(4,4,1)
    fTlocalx = floor(point[1])
    fTlocaly = floor(point[2])
    xBas0 = Int(fTlocalx)
    yBas0 = Int(fTlocaly)
    xBas1=xBas0+1
    yBas1=yBas0+1

    # Linear interpolation constants (percentages)
    xCom=point[1]-fTlocalx
    yCom=point[2]-fTlocaly
    xComi=(1-xCom)
    yComi=(1-yCom)
    perc[1]=xComi * yComi
    perc[2]=xComi * yCom
    perc[3]=xCom * yComi
    perc[4]=xCom * yCom

    # Stick to boundary
    if xBas0<1
        xBas0=1
        if xBas1<1
            xBas1=1
        end
    end

    if yBas0<1
        yBas0=1
        if yBas1<1
            yBas1=1
        end
    end
    if xBas1>Isize[1]
        xBas1=Isize[1]
        if xBas0>Isize[1]
            xBas0=Isize[1]
        end
    end

    if yBas1>Isize[2]
        yBas1=Isize[2]
        if yBas0>Isize[2]
            yBas0=Isize[2]
        end
    end

    # Get the neighbour intensities
    index[1] = LinearIndices(size(I))[xBas0,yBas0,1] #TODO Why is this dim 3,3,2?
    index[2] = LinearIndices(size(I))[xBas0,yBas1,1]
    index[3] = LinearIndices(size(I))[xBas1,yBas0,1]
    index[4] = LinearIndices(size(I))[xBas1,yBas1,1]

   f=Isize[1]*Isize[2]

    # the interpolated color
    Ireturn[1]=I[index[1]]*perc[1]+I[index[2]]*perc[2]+I[index[3]]*perc[3]+I[index[4]]*perc[4]
    Ireturn[2]=I[index[1]+f]*perc[1]+I[index[2]+f]*perc[2]+I[index[3]+f]*perc[3]+I[index[4]+f]*perc[4]

    return Ireturn
end

function rk4(startPoint::AbstractArray{T}, GradientVolume::AbstractArray{T}, Stepsize::Number) where T

    # Perform the RK4 raytracing step

   gradientArray = GradientVolume # TODO Test if this works (without copy)
    gradientArraySize = [size(gradientArray)...]
    nextPoint = startPoint

    startPoint1 = startPoint
    nextPoint =  RK4STEP_2D(gradientArray, startPoint1, ones(T,2,), Stepsize)
   return nextPoint
end


function RK4STEP_2D(gradientArray::AbstractArray{T}, startPoint::AbstractArray{T}, nextPoint,stepSize) where T

    # Perform one step of the RK4 algorithm
    k1 = ones(T,2,)
    k2 = ones(T,2,)
    k3 = ones(T,2,)
    k4 = ones(T,2,)
    tempPoint = ones(T,2,)

    gradientArraySize = [size(gradientArray)...]

    # Calculate k1
    k1 = interpgrad2d(k1, gradientArray, gradientArraySize, startPoint)
    tempnorm=max(norm(k1),1e-6)
    k1[1] = k1[1]*stepSize/tempnorm
    k1[2] = k1[2]*stepSize/tempnorm

    tempPoint[1]=startPoint[1] - k1[1]*0.5
    tempPoint[2]=startPoint[2] - k1[2]*0.5

    # Check if still inside the domain
    if !checkBounds2d([tempPoint...], gradientArraySize)
        return ones(2,)
    end

    # Calculate k2
    k2 = interpgrad2d(k2, gradientArray, gradientArraySize, tempPoint)
    tempnorm=max(norm(k2),1e-6)
    k2[1] = k2[1]*stepSize/tempnorm
    k2[2] = k2[2]*stepSize/tempnorm

    tempPoint[1]=startPoint[1] - k2[1]*0.5
    tempPoint[2]=startPoint[2] - k2[2]*0.5

    # Check the if are still inside the domain
    if !checkBounds2d([tempPoint...], gradientArraySize)
        return ones(2,)
    end

    # Calculate k3
    k3 = interpgrad2d(k3, gradientArray, gradientArraySize, tempPoint)
    tempnorm=max(norm(k3),1e-6)
    k3[1] = k3[1]*stepSize/tempnorm
    k3[2] = k3[2]*stepSize/tempnorm

    tempPoint[1]=startPoint[1] - k3[1]
    tempPoint[2]=startPoint[2] - k3[2]

    # Check the if are still inside the domain
    if !checkBounds2d([tempPoint...], gradientArraySize)
        return ones(2,)
    end

    # Calculate k4
    k4 = interpgrad2d(k4, gradientArray, gradientArraySize, tempPoint)
    tempnorm=max(norm(k4),1e-6)
    k4[1] = k4[1]*stepSize/tempnorm
    k4[2] = k4[2]*stepSize/tempnorm

    # Calculate final point
    nextPoint[1] = startPoint[1] - (k1[1] + k2[1]*2.0 + k3[1]*2.0 + k4[1])/6.0
    nextPoint[2] = startPoint[2] - (k1[2] + k2[2]*2.0 + k3[2]*2.0 + k4[2])/6.0

    # Set step to step size
    # /*
    # D[0]=(nextPoint[0]-startPoint[0]);
    # D[1]=(nextPoint[1]-startPoint[1]);
    # dl=stepSize/(sqrt(D[0]*D[0]+D[1]*D[1])+1e-15);
    # D[0]*=dl; D[1]*=dl;
    # nextPoint[0]=startPoint[0]+D[0];
    # nextPoint[1]=startPoint[1]+D[1];
    #*/

    # Check the if are still inside the domain
    if !checkBounds2d([nextPoint...], gradientArraySize)
        return ones(2,)
    end

    return nextPoint
end
