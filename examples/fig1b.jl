using CurveProximityQueries
using IntervalArithmetic
using Plots

pgfplots()

include("SymbolicCurves.jl")
using .SymbolicCurves

function heartcurve(a::Real, θ, xcenter, ycenter)
    x = :(($a*0.02*16)*(sin(t))^3)
    y = :(($a*0.02)*(13*cos(t) - 5*cos(t*2) - 2*cos(t*3) - cos(t*4)))
    SymCurve(x, y, 0..2π, θ, xcenter, ycenter)
end

function circinvolute(a::Real, θ, xcenter, ycenter)
    x = :($a*(cos(t) + t*sin(t)))
    y = :($a*(sin(t) - t*cos(t)))
    SymCurve(x, y, 0..10, θ, xcenter, ycenter)
end

a = heartcurve(1.0, 0.0, 0.0, 0.0)
b = circinvolute(0.035, -π, -1.0, 0.0)
pts = closest_points(a,b)
plot(a,
     linewidth=2.0,
     linecolor=:black,
    )
plot!(b,
      linewidth=2.0,
      linecolor=:black,
     )
plot!(pts,
      linestyle=:dash,
      markersize=8.0,
      linewidth=2.0,
     )
#  savefig("out/fig1b.pdf")
