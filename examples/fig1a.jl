using ConvexBodyProximityQueries
using CurveProximityQueries
using Plots
import Random: seed!

pgfplots()

seed!(22)
poly = randpoly([-0.5, 0.4]; n=20, scale=0.5)

curve = Bernstein([[0.944223, 0.177264], [0.46193, 0.105887], [0.791169, 0.172386], [0.256546, 0.463735], [0.70966, 0.3173], [0.664713, 0.830956], [0.952688, 0.736242], [0.00763912, 0.97189], [0.270595, 0.827442], [0.16636, 0.996615], [0.381073, 0.543077], [0.225381, 0.754049], [0.570043, 0.638543], [0.254696, 0.988563]])

pts = closest_points(poly, curve)

plot(poly,
     fillalpha=1.0,
     fillcolor=Gray(0.7),
    )
plot!(curve,
      linewidth=2.0,
      linecolor=:black,
     )
plot!(pts,
      linestyle=:dash,
      markersize=8.0,
      linewidth=2.0,
     )
savefig("out/fig1a.pdf")
