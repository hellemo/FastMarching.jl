using FastMarching
@static if VERSION >= v"0.7-"
  using Test
  using LinearAlgebra
else
  using Base.Test
end

function eye(n)
  fill(0,(n,n))+I
end


function testeye(n::Integer,endpoint=[1.,1.],stepsize=0.1)
  speedmap = eye(n) .* 1000 .+ 0.001
  source = [float(size(speedmap,1)), float(size(speedmap,2))]
  distancemap = FastMarching.msfm(speedmap,source,true,true)
  endpoint = [1.,1.]
  shortestline = FastMarching.shortestpath(distancemap,endpoint,source,0.1)
end

endpoint = [1.,1.]
@testset "testeye" begin
@testset "size: $t stepsize: $stepsize"  for t in 2:2:10, stepsize in [0.1,0.01]
  shortestline = testeye(t,endpoint,stepsize)
  @test sum(shortestline[:,2]-shortestline[:,1]) == 0
  @test sqrt((shortestline[1,1]-1.)^2-(shortestline[1,2]-1.)^2) <= stepsize
end
end

include(joinpath(@__DIR__, "multiplestarts.jl"))
