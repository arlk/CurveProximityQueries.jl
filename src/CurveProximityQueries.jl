module CurveProximityQueries

using LinearAlgebra
using StaticArrays
using IntervalArithmetic
using DataStructures
using Random: shuffle!
using CurveProximityQueries
import ConvexBodyProximityQueries: closest_points, minimum_distance, collision_detection, tolerance_verification

export closest_points, minimum_distance, collision_detection, tolerance_verification

include("parametric.jl")
include("bounds.jl")
include("bnb.jl")
include("obstacles.jl")
include("plot.jl")

end # module
