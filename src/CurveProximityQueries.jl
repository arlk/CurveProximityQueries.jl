module CurveProximityQueries

export closest_points, minimum_distance, collision_detection, tolerance_verification

using LinearAlgebra
using StaticArrays
using IntervalArithmetic
using DataStructures
using Random: shuffle!
import ConvexBodyProximityQueries

include("parametric.jl")
include("bounds.jl")
include("bnb.jl")
include("obstacles.jl")
include("plot.jl")

end # module
