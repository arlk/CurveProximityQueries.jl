export Curve
export Bernstein
export differentiate, integrate, arclength

abstract type Curve{D, T} end

include("curves/bernstein.jl")
