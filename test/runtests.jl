using FastMarching
using Base.Test
using LinearAlgebra

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

function runtests()
  const endpoint = [1.,1.]
  for t in 2:2:10, stepsize in [0.1,0.01]
    println("size: $(lpad(t,2," ")) stepsize: $(rpad(string(stepsize),4,"0"))")
    shortestline = testeye(t,endpoint,stepsize)
    @test sum(shortestline[:,2]-shortestline[:,1]) == 0
    @test sqrt((shortestline[1,1]-1.)^2-(shortestline[1,2]-1.)^2) <= stepsize
  end

  try
    rm("multiplestarts.png")
  end
  include(joinpath(Pkg.dir("FastMarching"),"examples/multiplestart.jl"))
  @test isfile("multiplestart.png")


end

runtests()
