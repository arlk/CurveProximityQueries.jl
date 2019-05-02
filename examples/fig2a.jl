using CurveProximityQueries
using Plots
import Random: seed!

pgfplots()
seed!(19)

a = rand(Bernstein{2,12})
cp = a.f.control_points .- 0.5
a = Bernstein(cp..., 2*cp[1]-cp[2], cp[1])

b = Bernstein([[0.95, -0.2], [0.85, -0.2], [0.85, -0.1], [0.9125, -0.1], [0.975, -0.1], [1.0375, -0.1], [1.1, -0.1], [1.1, -0.1875], [1.1, -0.275], [1.1, -0.3625], [1.1, -0.45], [0.975, -0.45], [0.85, -0.45], [0.725, -0.45], [0.6, -0.45], [0.6, -0.3375], [0.6, -0.225], [0.6, -0.1125], [0.6, 0.0], [0.675, 0.0], [0.75, 0.0], [0.825, 0.0], [0.9, 0.0], [0.9, 0.1], [0.8, 0.1], [0.75, 0.1], [0.7, 0.1], [0.7, 0.05], [0.75, 0.0], [0.8, -0.05], [0.9, -0.05]])

pts = closest_points(a, b)

plot(a,
     linecolor=:black,
     linewidth=1.0,
     fill=(0, Gray(0.7)),
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
annotate!(0.07, -0.2, text("\$\\mathcal{B}\$", 25))
annotate!(0.8, -0.2, text("\$\\Psi_{\\mathcal{I}}\$", 25))
annotate!(pts[1][1]+0.05, pts[1][2]+0.05, text("\$b^*\$", 22))
annotate!(pts[2][1]-0.10, pts[2][2]+0.05, text("\$\\psi(t^*)\$", 22))
savefig("out/fig2a.pdf")
