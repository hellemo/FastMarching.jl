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

function checkBounds2d( point::Array{Float64,1}, Isize::Array{Int,1})
    if (point[1]<1)||(point[2]<1)||(point[1]>(Isize[1]))||(point[2]>(Isize[2]))
        return false
    end
    return true
end

# TODO check dimensions of point, seems to be bogus
function checkBounds3d( point::Array{Float64,3}, Isize::Array{Int,1})
    if (point[1,1,1]<1)||(point[2,1,1]<1)||(point[3,1,1]<1)||(point[1,1,1]>(Isize[1]))||(point[2,1,1]>(Isize[2]))||(point[3,1,1]>(Isize[3]))
         return false
    end
    return true
end

# TODO: Replace with (wrapped) built in norm?
function norm2(a::Array{Float64,1})
    return sqrt(a[1]*a[1]+a[2]*a[2])
end

# TODO: Replace with (wrapped) built in norm?
function norm3(a::Array{Float64,1})
    return sqrt(a[1]*a[1]+a[2]*a[2]+a[3]*a[3])
end


#function interpgrad2d(k1, gradientArray, gradientArraySize, startPoint,tmp)
#    tmp = e1(copy(startPoint),gradientArray,0.5)
#    k1[1] = tmp[1]
#    k1[2] = tmp[2]
#end

function interpgrad2d(Ireturn::Array{Float64,1},I::Array{Float64,3}, Isize::Array{Int64,1}, point::Array{Float64,1})
    #  Linear interpolation variables

    perc = zeros(4,)
    #index = [1::Int64 for i = 1:4]
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
# TODO: Port 3D
# __inline void interpgrad3d(double *Ireturn, double *I, int *Isize, double *point) {
#     /*  Linear interpolation variables */
#     int xBas0, xBas1, yBas0, yBas1, zBas0, zBas1;
#     double perc[8];
#     double xCom, yCom, zCom;
#     double xComi, yComi, zComi;
#     double fTlocalx, fTlocaly, fTlocalz;
#     int f0, f1;
#     int index[8];
#     double temp;

#     fTlocalx = floor(point[0]); fTlocaly = floor(point[1]); fTlocalz = floor(point[2]);
#     xBas0=(int) fTlocalx; yBas0=(int) fTlocaly; zBas0=(int) fTlocalz;
#     xBas1=xBas0+1; yBas1=yBas0+1; zBas1=zBas0+1;

#     /* Linear interpolation constants (percentages) */
#     xCom=point[0]-fTlocalx;  yCom=point[1]-fTlocaly;   zCom=point[2]-fTlocalz;
#     xComi=(1-xCom); yComi=(1-yCom); zComi=(1-zCom);
#     perc[0]=xComi * yComi; perc[1]=perc[0] * zCom; perc[0]=perc[0] * zComi;
#     perc[2]=xComi * yCom;  perc[3]=perc[2] * zCom; perc[2]=perc[2] * zComi;
#     perc[4]=xCom * yComi;  perc[5]=perc[4] * zCom; perc[4]=perc[4] * zComi;
#     perc[6]=xCom * yCom;   perc[7]=perc[6] * zCom; perc[6]=perc[6] * zComi;

#     /* Stick to boundary */
#     if(xBas0<0) { xBas0=0; if(xBas1<0) { xBas1=0; }}
#     if(yBas0<0) { yBas0=0; if(yBas1<0) { yBas1=0; }}
#     if(zBas0<0) { zBas0=0; if(zBas1<0) { zBas1=0; }}
#     if(xBas1>(Isize[0]-1)) { xBas1=Isize[0]-1; if(xBas0>(Isize[0]-1)) { xBas0=Isize[0]-1; }}
#     if(yBas1>(Isize[1]-1)) { yBas1=Isize[1]-1; if(yBas0>(Isize[1]-1)) { yBas0=Isize[1]-1; }}
#     if(zBas1>(Isize[2]-1)) { zBas1=Isize[2]-1; if(zBas0>(Isize[2]-1)) { zBas0=Isize[2]-1; }}

#    /*Get the neighbour intensities */
#     index[0]=mindex3(xBas0, yBas0, zBas0, Isize[0], Isize[1]);
#     index[1]=mindex3(xBas0, yBas0, zBas1, Isize[0], Isize[1]);
#     index[2]=mindex3(xBas0, yBas1, zBas0, Isize[0], Isize[1]);
#     index[3]=mindex3(xBas0, yBas1, zBas1, Isize[0], Isize[1]);
#     index[4]=mindex3(xBas1, yBas0, zBas0, Isize[0], Isize[1]);
#     index[5]=mindex3(xBas1, yBas0, zBas1, Isize[0], Isize[1]);
#     index[6]=mindex3(xBas1, yBas1, zBas0, Isize[0], Isize[1]);
#     index[7]=mindex3(xBas1, yBas1, zBas1, Isize[0], Isize[1]);
#     f0=Isize[0]*Isize[1]*Isize[2];
#     f1=f0+f0;

#    /*the interpolated color */
#     temp=I[index[0]]*perc[0]+I[index[1]]*perc[1]+I[index[2]]*perc[2]+I[index[3]]*perc[3];
#     Ireturn[0]=temp+I[index[4]]*perc[4]+I[index[5]]*perc[5]+I[index[6]]*perc[6]+I[index[7]]*perc[7];
#     temp=I[index[0]+f0]*perc[0]+I[index[1]+f0]*perc[1]+I[index[2]+f0]*perc[2]+I[index[3]+f0]*perc[3];
#     Ireturn[1]=temp+I[index[4]+f0]*perc[4]+I[index[5]+f0]*perc[5]+I[index[6]+f0]*perc[6]+I[index[7]+f0]*perc[7];
#     temp=I[index[0]+f1]*perc[0]+I[index[1]+f1]*perc[1]+I[index[2]+f1]*perc[2]+I[index[3]+f1]*perc[3];
#     Ireturn[2]=temp+I[index[4]+f1]*perc[4]+I[index[5]+f1]*perc[5]+I[index[6]+f1]*perc[6]+I[index[7]+f1]*perc[7];
# }
function rk4(startPoint::Array{Float64,1}, GradientVolume::Array{Float64,3}, Stepsize::Float64)

    # Perform the RK4 raytracing step

    #gradientArraySize[0]=gradientArraySizeC[0];
    #    gradientArraySize[1]=gradientArraySizeC[1];



   #gradientArray = copy(GradientVolume)
   gradientArray = GradientVolume # TODO Test if this works (without copy)
    gradientArraySize = [size(gradientArray)...]
    nextPoint = startPoint

    startPoint1 = startPoint
#    startPoint1[1]=startPoint[1]-1.0
#    startPoint1[2]=startPoint[2]-1.0
#    if RK4STEP_2D(gradientArray, gradientArraySize, startPoint1, nextPoint, Stepsize)
    nextPoint =  RK4STEP_2D(gradientArray, startPoint1, ones(2,), Stepsize)
 #       nextPoint[1]=nextPoint[1]+1.0
 #       nextPoint[2]=nextPoint[2]+1.0
#    else
#        nextPoint[1]=1.
#        nextPoint[2]=1.
 #   end
    # TODO: 3D
    # else if(PointLength==3) {
    #     gradientArraySize[0]=gradientArraySizeC[0];
    #     gradientArraySize[1]=gradientArraySizeC[1];
    #     gradientArraySize[2]=gradientArraySizeC[2];
    #     startPoint1[0]=startPoint[0]-1.0;
    #     startPoint1[1]=startPoint[1]-1.0;
    #     startPoint1[2]=startPoint[2]-1.0;
    #     if(RK4STEP_3D(gradientArray, gradientArraySize, startPoint1, nextPoint, stepSize)) {
    #         nextPoint[0]=nextPoint[0]+1.0;
    #         nextPoint[1]=nextPoint[1]+1.0;
    #         nextPoint[2]=nextPoint[2]+1.0;
    #     }
    #     else {
    #         nextPoint[0]=0; nextPoint[1]=0; nextPoint[2]=0;
    #     }

    # }
    # else {
    #     mexErrMsgTxt("Starting Point must be 2D or 3D");
    # }

   return nextPoint
end


function RK4STEP_2D(gradientArray::Array{Float64,3}, startPoint::Array{Float64,1}, nextPoint::Array{Float64,1},stepSize::Float64)

    # Perform one step of the RK4 algorithm
    k1 = ones(2,)
    k2 = ones(2,)
    k3 = ones(2,)
    k4 = ones(2,)
    tempPoint = ones(2,)

    gradientArraySize = [size(gradientArray)...]

    # Calculate k1
    k1 = interpgrad2d(k1, gradientArray, gradientArraySize, startPoint)
    tempnorm=max(norm2(k1),1e-6)
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
    tempnorm=max(norm2(k2),1e-6)
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
    tempnorm=max(norm2(k3),1e-6)
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
    tempnorm=max(norm2(k4),1e-6)
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
