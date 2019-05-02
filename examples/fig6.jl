using ConvexBodyProximityQueries
using CurveProximityQueries
using Plots
import Random: seed!

pgfplots()

include("EllipseBounds.jl")
using .EllipseBounds

function drawellipses!(a, b::Curve, maxnum; atol=1e-8)
    leaves, glb, gub, soln = CurveProximityQueries._init(a, b)
    i = 0
    while (gub - glb) > atol && i < maxnum
        glb, gub, soln = CurveProximityQueries._loop(a, b, leaves, glb, gub)
        i += 1
    end
    for interval in leaves
        β, lb = interval
        if diam(β) > 1e-3
            bell = Ellipse(b, β)
            plot!(bell,
                 linecolor=RGBA(102/255,130/225,223/255,0.5),
                 fillcolor=RGBA(102/255,130/225,223/255,0.5),
                 )
        end
    end
end

seed!(22)
poly = randpoly([-0.5, 0.4]; n=20, scale=0.5)

curve = Bernstein([[0.944223, 0.177264], [0.46193, 0.105887], [0.791169, 0.172386], [0.256546, 0.463735], [0.70966, 0.3173], [0.664713, 0.830956], [0.952688, 0.736242], [0.00763912, 0.97189], [0.270595, 0.827442], [0.16636, 0.996615], [0.381073, 0.543077], [0.225381, 0.754049], [0.570043, 0.638543], [0.254696, 0.988563]])

pts = closest_points(poly, curve)

cases = [nothing, 0, 1, 3, 7, 100]
alph = 'a':'z'
for (i, c) in enumerate(cases)
    plot(poly,
         fillalpha=1.0,
         fillcolor=Gray(0.7),
        )
    isnothing(c) || drawellipses!(poly, curve, c)
    plot!(curve,
          linewidth=2.0,
          linecolor=:black,
         )
    savefig("out/fig6$(alph[i]).pdf")
end
plot!(pts,
      linestyle=:dash,
      markersize=8.0,
      linewidth=2.0,
     )
savefig("out/fig6$(alph[length(cases)]).pdf")
