# FastMarching

This is a (partial) port of  [Accurate Fast Marching](https://se.mathworks.com/matlabcentral/fileexchange/24531-accurate-fast-marching?ue) by Dirk-Jan Kroon

## Example

Multiple starting points:
```julia
using FastMarching, Gadfly
tsize = 300
npoints = 10
SpeedImage = ones(tsize,tsize)
SourcePoint = rand(2,npoints).*tsize.+1

t1 = FastMarching.msfm(SpeedImage, SourcePoint, true, true)

draw(SVG("examples/multistart.svg"),
    Gadfly.plot(z=t1, Geom.contour))

```
![](examples/multistart.svg)


Originally ported to Julia 0.2, updated for 1.0 2018.

[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
