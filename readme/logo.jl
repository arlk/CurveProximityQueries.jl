using CurveProximityQueries
using ConvexProximityQueries
using StaticArrays
using Plots
using Random: seed!

function movecp(t; k=3)
   x1 = 3*cos(k*t)*cos(t)
   y1 = cos(k*t)*sin(t)
   x2 = 3*(sin(t)+1.0)
   y2 = 3*cos(t)
   SVector{2}(x1, y1), SVector{2}(x2, y2)
end

anim = @animate for i=1:360
    t = i/180*π
    b = [SVector{2}(0.0, 0.0),
         movecp(t)...,
         SVector{2}(3.0, 0.0)]
    B = Bernstein(b...)

    seed!(1)
    A = randpoly([1.5, 0.0], t; scale=1.0)

    plot(A)
    plot!(xlims = (-1.5,4.5), ylims = (-2.0,2.0))

    @show t
    if collision_detection(A, B)
        plot!(B, linecolor = RGB(213/255,99/255,92/255), linewidth = 4)
    else
        AB = closest_points(A, B)
        plot!(AB, linewidth = 3)
        plot!(B, linecolor = RGB(106/255,92/225,213/255), linewidth = 4)
    end
end

gif(anim, "logo1.gif", fps = 30)

anim = @animate for i=1:360
    t = i/180*π
    b = [SVector{2}(0.0, 0.0),
         movecp(t)...,
         SVector{2}(3.0, 0.0)]
    B = Bernstein(b...)

    a = [SVector{2}(0.0, 1.5),
         movecp(t; k=2)...,
         SVector{2}(3.0, 1.5)]
    A = Bernstein(a...)

    plot(xlims = (-1.5,4.5), ylims = (-2.0,2.0))

    @show t
    if collision_detection(A, B; atol=1e-7)
        plot!(B, linecolor = RGB(213/255,99/255,92/255), linewidth = 4)
        plot!(A, linecolor = RGB(213/255,99/255,92/255), linewidth = 4)
    else
        AB = closest_points(A, B; atol=1e-7)
        plot!(AB, linewidth = 3)
        plot!(B, linecolor = RGB(106/255,92/225,213/255), linewidth = 4)
        plot!(A, linecolor = RGB(120/255,213/255,104/255), linewidth = 4)
    end
end

gif(anim, "logo2.gif", fps = 30)
