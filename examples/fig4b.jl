using CurveProximityQueries
using ConvexBodyProximityQueries
using IntervalArithmetic
using Plots
using StaticArrays

pgfplots()

include("EllipseBounds.jl")
using .EllipseBounds

curve = Bernstein([[-0.6, 0.0], [-0.6, 0.0], [0.1, -0.4], [0.4, 1.6], [-0.2, -1.6], [-0.1, 0.4], [0.8, 0.0], [0.8, 0.0]])
ell = Ellipse(curve, 0.2..0.8)

curve2 = Bernstein([[0.6, 0.6], [0.4, 1.0], [0.2, 0.2], [0.0, 0.6]])
ell2 = Ellipse(curve2, 0.2..0.9)

pts = closest_points(ell, ell2, @SVector(ones(2)))

plot(ell,
    linecolor=RGBA(102/255,130/225,223/255,0.5),
    fillcolor=RGBA(102/255,130/225,223/255,0.5),
    )
plot!(ell2,
     linecolor=RGBA(102/255,130/225,223/255,0.5),
     fillcolor=RGBA(102/255,130/225,223/255,0.5),
     )
plot!(curve, 0.0..0.2,
      linewidth=2.0,
      linestyle=:dash,
      linecolor=Gray(0.5),
     )
plot!(curve, 0.2..0.8,
      linewidth=2.0,
      linecolor=:black,
     )
plot!(curve, 0.8..1.0,
      linewidth=2.0,
      linestyle=:dash,
      linecolor=Gray(0.5),
     )
scatter!(curve(0.2), markersize=4, markercolor=:black)
scatter!(curve(0.8), markersize=4, markercolor=:black)

plot!(curve2, 0.0..0.2,
      linewidth=2.0,
      linestyle=:dash,
      linecolor=Gray(0.5),
     )
plot!(curve2, 0.2..0.9,
      linewidth=2.0,
      linecolor=:black,
     )
plot!(curve2, 0.9..1.0,
      linewidth=2.0,
      linestyle=:dash,
      linecolor=Gray(0.5),
     )
scatter!(curve2(0.2), markersize=4, markercolor=:black)
scatter!(curve2(0.9), markersize=4, markercolor=:black)

plot!(pts,
      linestyle=:dash,
      markersize=8.0,
      linewidth=2.0,
     )
savefig("out/fig4b.pdf")
