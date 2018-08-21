using FastMarching
using Images
using FileIO

function maze()
  Float64.(channelview(img))
  img = load(joinpath(Pkg.dir("FastMarching"),"examples/images/maze.png"))
end

maze()
