using ConvexBodyProximityQueries
using CurveProximityQueries
using Plots
import Random: seed!

pgfplots()

include("EllipseBounds.jl")
using .EllipseBounds

function geterr(a, b::Curve; atol=1e-10)
    err = []
    leaves, glb, gub, soln = CurveProximityQueries._init(a, b)
    i = 0
    push!(err, gub-glb)
    while (gub - glb) > atol
        glb, gub, soln = CurveProximityQueries._loop(a, b, leaves, glb, gub)
        push!(err, gub-glb)
        i += 1
    end
    print("\nTotal Iter: $i\n")
    return err
end

function getLogTicks(x)
    min = ceil(log10(minimum(x)))
    max = ceil(log10(maximum(x)))
    major = 10 .^collect(min:max)
    majorText = []
    for i = min:max
        if i == 0
            push!(majorText, "1")
        elseif i%2 == 0
            push!(majorText, "10^{$(round(Int64,i))")
        else
            push!(majorText, "")
        end
    end
    minor = [j*10^i for i=(min-1):(max+1) for j=2:9]
    minor = minor[findall(minimum(x) .<= minor .<= maximum(x))]
    ([major; minor] , [majorText; fill("", length(minor))])
end

seed!(22)
poly = randpoly([-0.5, 0.4]; n=20, scale=0.5)

curve = Bernstein([[0.944223, 0.177264], [0.46193, 0.105887], [0.791169, 0.172386], [0.256546, 0.463735], [0.70966, 0.3173], [0.664713, 0.830956], [0.952688, 0.736242], [0.00763912, 0.97189], [0.270595, 0.827442], [0.16636, 0.996615], [0.381073, 0.543077], [0.225381, 0.754049], [0.570043, 0.638543], [0.254696, 0.988563]])

err = geterr(poly, curve)

tic = getLogTicks(err)
plot(err, linetype=:steppost, linecolor=:black, linewidth=1.5, yscale=:log10, yticks=tic)
plot!(legend=false, xlabel="\$\\mathrm{card}(\\mathcal{L})\$", ylabel="\$\\overline{d}- \\underline{d}\$", guidefontsize=18)
plot!(gridlinewidth = 0.5, gridstyle = :solid, gridalpha = 0.2)
savefig("out/fig7.pdf")
