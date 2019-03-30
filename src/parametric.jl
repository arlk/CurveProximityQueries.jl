export ParametricCurve, RectifiableCurve
export Bernstein, @bernstein
export RectifiableBernstein, @rectifiablebernstein
export differentiate, integrate, arclength

abstract type ParametricCurve{D, T} end
abstract type RectifiableCurve{D, T} end

include("curves/bernstein.jl")
