module CurveProximityQueries

using LinearAlgebra
using StaticArrays
using IntervalArithmetic
using DataStructures
import ConvexBodyProximityQueries: closest_points, minimum_distance, collision_detection, tolerance_verification

export closest_points, minimum_distance, collision_detection, tolerance_verification

include("parametric.jl")
include("bounds.jl")
include("bnb.jl")
include("plot.jl")

end # module
