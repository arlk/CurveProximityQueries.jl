# CurveProximityQueries

[![Build Status](https://travis-ci.org/arlk/CurveProximityQueries.jl.svg?branch=master)](https://travis-ci.org/arlk/CurveProximityQueries.jl) [![codecov](https://codecov.io/gh/arlk/ConvexBodyProximityQueries.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/arlk/CurveProximityQueries.jl)

Curve - Convex Polygon        |  Curve - Curve
:-------------------------:|:-------------------------:
![](https://github.com/arlk/CurveProximityQueries.jl/raw/master/readme/logo1.gif)  |  ![](https://github.com/arlk/CurveProximityQueries.jl/raw/master/readme/logo2.gif)

CurveProximityQueries implements methods to compute the
 - Closest Points
 - Minimum Distance
 - Tolerance Verification
 - Collision Detection

between an absolutely continuous parametric curve and another object, or between two curves. If you find this package/work useful in your research please cite our [paper](http://www.roboticsproceedings.org/rss15/p42.pdf):
```
@INPROCEEDINGS{Hovakimyan-RSS-19, 
    AUTHOR    = {Arun Lakshmanan AND Andrew Patterson AND Venanzio Cichella AND Naira Hovakimyan}, 
    TITLE     = {Proximity Queries for Absolutely Continuous Parametric Curves}, 
    BOOKTITLE = {Proceedings of Robotics: Science and Systems}, 
    YEAR      = {2019}, 
    ADDRESS   = {Freiburg im Breisgau, Germany}, 
    MONTH     = {June}, 
    DOI       = {10.15607/RSS.2019.XV.042} 
}
```

## Usage

### Installation

```julia
julia> ] add CurveProximityQueries
```

### Curve Types

Methods are available to create and manipulate Bernstein polynomials. A Bernstein polynomial with uniformly randomly sampled control points can be created with:
```julia
julia> rand(Bernstein{3, 6})
a 5th order Bernstein polynomial with control points at:
([0.345747, 0.149338, 0.609345], [0.86539, 0.736102, 0.31424], [0.20149, 0.167441, 0.662126], [0.975667, 0.468056, 0.505278], [0.371533, 0.477992, 0.83668], [0.322821, 0.973494, 0.93129])
with an arclength of 1.463398157000094
```
which constructs a 3D 5th order Bernstein polynomial. Control points can be directly fed to the constructor as well:
```julia
julia> cpts = [[0.0, 0.0], [1.0, 2.0], [2.0, 0.0], [3.0, 0.0]];
julia> c = Bernstein(cpts)
a 3rd order Bernstein polynomial with control points at:
([0.0, 0.0], [1.0, 2.0], [2.0, 0.0], [3.0, 0.0])
with an arclength of 3.714835124201342
```
Generally Bernstein polynomials are evaluated between $[0, 1]$, but custom limits can be provided using `Interval` from the `IntervalArithmetic` package:
```julia
julia> c = Bernstein(cpts, limits=Interval(0.5, 2.5))
a 3rd order Bernstein polynomial with control points at:
([0.0, 0.0], [1.0, 2.0], [2.0, 0.0], [3.0, 0.0])
with an arclength of 3.714835124201342
```
The `Bernstein` types are functors that evaluate the polynomial at some value.
```julia
julia> c(1.5)
2-element SArray{Tuple{2},Float64,1,2}:
 1.5
 0.75
```
Note, that the arclength is not a true arclength but an upper bound that is used by the proximity methods. This upper bound can be evaluated at an interval or from the beginning of the limit using:
```julia
julia> arclength(c, Interval(0.9, 1.4))
0.7839413641976036
julia> arclength(c, 1.0) # evaluates from 0.5 to 1.0
1.1826380733343569
```
Several algebraic operations over Bernstein polynomials are available for convenience: addition, subtraction, product with reals, differentiation.

### Obstacle Types

Several obstacle types are provided with convenience macros from [ConvexBodyProximityQueries.jl](https://github.com/arlk/ConvexBodyProximityQueries.jl):
```julia
julia> @point rand(3)
ConvexPolygon{3,1,Float64}(SArray{Tuple{3},Float64,1,3}[[0.135678, 0.840508, 0.140532]])
julia> @line [0.,1.,1.], [1.,2.,1.] # point A, point B
ConvexPolygon{3,2,Float64}(SArray{Tuple{3},Float64,1,3}[[0.0, 1.0, 1.0], [1.0, 2.0, 1.0]])
julia> @rect [0.,0.], rand(2) # center, widths
ConvexPolygon{2,4,Float64}(SArray{Tuple{2},Float64,1,2}[[0.395191, 0.174093], [-0.395191, 0.174093], [-0.395191, -0.174093], [0.395191, -0.174093]])
julia> @square ones(3), 1.0 # center, width
ConvexPolygon{3,8,Float64}(SArray{Tuple{3},Float64,1,3}[[1.5, 1.5, 1.5], [0.5, 1.5, 1.5], [0.5, 0.5, 1.5], [1.5, 0.5, 1.5], [1.5, 1.5, 0.5], [0.5, 1.5, 0.5], [0.5, 0.5, 0.5], [1.5, 0.5, 0.5]])
```
Random convex polygons can be constructed for 2D:
```julia
julia> obs = randpoly([1., 2.], 0.5; scale=1.0, n=20) # center, rotation; scale, number of vertices
ConvexPolygon{2,20,Float64}(SArray{Tuple{2},Float64,1,2}[[0.642686, 2.36248], [0.622121, 2.34973], [0.42866, 2.06399], [0.412454, 2.0344], [0.454968, 1.98069], [0.499506, 1.92797], [0.599317, 1.82251], [0.62982, 1.79366], [0.659987, 1.76526], [0.733777, 1.71118], [0.87861, 1.63702], [1.07313, 1.54129], [1.46142, 1.68951], [1.46817, 1.72673], [1.48588, 1.85669], [1.46772, 2.06245], [1.3987, 2.23026], [1.30631, 2.4218], [1.20662, 2.61294], [0.88346, 2.47282]])
```
### Proximity Queries
The formal definitions for each of the proximity queries can be found in the paper.

#### Minimum Distance
```julia
julia> minimum_distance(c, obs)
0.6542938586889715
```

#### Tolerance Verification
```julia
julia> tolerance_verification(c, obs, 0.5)
true
julia> tolerance_verification(c, obs, 1.0)
false
```

#### Collision Detection
```julia
julia> collision_detection(c, obs)
false
```

#### Closest Points
```julia
julia> pts = closest_points(c, obs)
([1.07313, 1.54129], [1.03968, 0.887853])
```

All the above tests can be performed when `obs` is replaced with another parametric curve.

### Plotting
Plotting recipes are provided for the curves. For example to plot the closest points between `obs` and `c`, one can simply:
```julia
julia> plot(c); plot!(obs); plot!(pts)
```
![](https://github.com/arlk/CurveProximityQueries.jl/raw/master/readme/example.png)

### Custom Types

Parametric curves can be user defined. For example, a monomial basis type can be created as a subtype of `Curve{D, T}` where `D` is the dimension of the curve and `T` is the eltype:
```julia
struct Monomial{D, N, T} <: Curve{D, T}
  coeff::SVector{N, SVector{D, T}}
  limits::Interval{T}
end
```
Methods to compute the evaluation (as a functor), and arclength upper bound has to be provided (see paper for details).

Similarly, custom types for obstacles can be created. If the obstacles are convex, then only a [support mapping](https://github.com/arlk/ConvexBodyProximityQueries.jl#usage) for the custom type is required. However, if the obstacle is not convex then methods to compute the `minimum_distance`, `tolerance_verification`, and `collision_detection` between a point in space and the object must be provided.
