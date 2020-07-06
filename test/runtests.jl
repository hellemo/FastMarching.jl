using FastMarching
using Test
using LinearAlgebra

@testset "Test FastMarching" begin

function testI(n::Integer,endpoint=[1.,1.],stepsize=0.1,T=Float64)
  speedmap = T.(fill(0,(n,n))+I .* 1000 .+ 0.1)
  source = T.([float(size(speedmap,1)), float(size(speedmap,2))])
  distancemap = FastMarching.msfm(speedmap,source,true,true)
  endpoint = T.([1.,1.])
  stepsize = T(stepsize)
  shortestline = FastMarching.shortestpath(distancemap,endpoint,source,stepsize)
end

endpoint = [1.,1.]
@testset "testeye" begin
@testset "size: $t stepsize: $stepsize"  for t in 2:2:10, stepsize in [0.1,0.01], T in [Float16,Float32,Float64]
  shortestline = testI(t,endpoint,stepsize,T)
  @test sum(shortestline[:,2]-shortestline[:,1]) == 0
  @test sqrt((shortestline[1,1]-1.)^2-(shortestline[1,2]-1.)^2) <= stepsize
  @test eltype(shortestline) == T
end
end

include(joinpath(@__DIR__, "multiplestarts.jl"))
end