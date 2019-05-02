using CurveProximityQueries
using IntervalArithmetic
using Plots

pgfplots()

include("EllipseBounds.jl")
using .EllipseBounds

curve = Bernstein([[-0.6, 0.0], [-0.6, 0.0], [0.1, -0.4], [0.4, 1.6], [-0.2, -1.6], [-0.1, 0.4], [0.8, 0.0], [0.8, 0.0]])
ell = Ellipse(curve, 0.2..0.8)

plot(ell,
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
annotate!(curve(0.5)[1]+0.03, curve(0.5)[2]+0.17, "\$\\mathcal{U}_{\\mathcal{Q}}\$")
annotate!(curve(0.2)[1], curve(0.2)[2]-0.07, "\$\\psi(\\alpha)\$")
annotate!(curve(0.8)[1], curve(0.8)[2]+0.07, "\$\\psi(\\beta)\$")
savefig("out/fig3.pdf")
