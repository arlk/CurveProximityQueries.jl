using ConcaveHull
using CurveProximityQueries
using LinearAlgebra
using Plots
import Random: seed!

pgfplots()
seed!(19)

@recipe function f(a::ConcaveHull.Hull)
    grid := false
    axis := false
    ticks := false
    legend := false
    aspect_ratio := :equal
    linewidth --> 1.5
    linecolor --> RGBA(220/255, 20/255, 60/255, 0.6) #Gray(0.7)
    linestyle --> :dash
    leg := false
    fill --> (0, RGBA(220/255, 20/255, 60/255, 0.1))
    @series begin
        pts = hcat(a.vertices...)
        x = vcat(pts[1, :], pts[1, 1])
        y = vcat(pts[2, :], pts[2, 1])
        x, y
    end
end

a = rand(Bernstein{2,12})
cp = a.f.control_points .- 0.5
a = Bernstein(cp..., 2*cp[1]-cp[2], cp[1])

b = Bernstein([[0.95, -0.2], [0.85, -0.2], [0.85, -0.1], [0.9125, -0.1], [0.975, -0.1], [1.0375, -0.1], [1.1, -0.1], [1.1, -0.1875], [1.1, -0.275], [1.1, -0.3625], [1.1, -0.45], [0.975, -0.45], [0.85, -0.45], [0.725, -0.45], [0.6, -0.45], [0.6, -0.3375], [0.6, -0.225], [0.6, -0.1125], [0.6, 0.0], [0.675, 0.0], [0.75, 0.0], [0.825, 0.0], [0.9, 0.0], [0.9, 0.1], [0.8, 0.1], [0.75, 0.1], [0.7, 0.1], [0.7, 0.05], [0.75, 0.0], [0.8, -0.05], [0.9, -0.05]])

points = []
da = CurveProximityQueries.differentiate(a.f)
for v = 0:0.002:1
    pt = a(v)
    mvec = [0. 1; -1 0]*normalize(da(v))
    push!(points, Array(pt+mvec*0.13))
end
ahull = concave_hull(points)

plot(ahull)
plot!(a,
      linecolor=:black,
      linewidth=1.0,
      fill=(0, Gray(0.7)),
     )
plot!([0.282932, 0.412878], [-0.195185, -0.198932], linecolor = Gray(0.2), markershape = :circle, markersize = 4, markercolor = :black, markerstrokecolor = :black)
plot!(b,
      linewidth=2.0,
      linecolor=:black,
     )
annotate!(0.07, -0.2, text("\$\\mathcal{B}\$", 25))
annotate!(0.8, -0.2, text("\$\\Psi_{\\mathcal{I}}\$", 25))
annotate!(0.348, -0.147, text("\$\\Delta\$", 22))
savefig("out/fig2b.pdf")
