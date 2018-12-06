using Gadfly


function testMultipleStarts(tsize::Int=300,npoints::Int=20)
    SpeedImage = ones(tsize,tsize)
    SourcePoint = rand(2,npoints).*tsize.+1
    FastMarching.msfm(SpeedImage,SourcePoint,true,true)
end  # function testMultipleStarts

t1 = testMultipleStarts()
draw(PNG("multiplestart.png"),plot(z=t1,Geom.contour))