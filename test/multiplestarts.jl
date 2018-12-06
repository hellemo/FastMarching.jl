using FastMarching
@static if VERSION >= v"0.7-"
  using Test
else
  using Base.Test
end

function testMultipleStarts(tsize::Int=300,npoints::Int=20)
  SpeedImage = ones(tsize,tsize)
  SourcePoint = rand(2,npoints).*tsize.+1

  t1 = FastMarching.msfm(SpeedImage, SourcePoint, true, true)

  t2 = FastMarching.msfm(SpeedImage, SourcePoint, true, true)
  @test norm(t1 - t2) < 1e-10 # Result is deterministic
end  # function testMultipleStarts

@testset "testMultipleStarts is deterministic" begin
@testset "npoints=$npoints" for npoints in [10, 15]
  testMultipleStarts(300, npoints)
end
end