export Curve
export Bernstein
export arclength

abstract type Curve{D, T} end

include("curves/bernstein.jl")
