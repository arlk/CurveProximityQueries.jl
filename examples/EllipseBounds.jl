module EllipseBounds
using CurveProximityQueries
using ConvexBodyProximityQueries
using IntervalArithmetic
using LinearAlgebra
using RecipesBase
using StaticArrays

export Ellipse

struct Ellipse{A, B, C}
    center::A
    rotate::C
    radii::B
end

import ConvexBodyProximityQueries: support
function support(el::Ellipse, dir::SVector)
    rotdir = el.rotate'*dir
    sp = el.radii .* normalize(el.radii .* rotdir)
    el.rotate*sp .+ el.center
end

function Ellipse(c::Curve, θ::Interval)
    focusA = c(θ.lo)
    focusB = c(θ.hi)
    center = (focusB + focusA)/2.0
    BA = focusB - focusA
    focaldistance = norm(BA)
    dir = BA/focaldistance
    rotate = rotationmatrix(dir)
    majradius = arclength(c, θ)/2
    minradius = majradius > focaldistance/2 ? sqrt(majradius^2 - (focaldistance/2)^2) : 0.0
    radii = SVector{2}(majradius, minradius)
    Ellipse(center, rotate, radii)
end

function rotationmatrix(dir::SVector{2})
    SMatrix{2,2}(dir[1],dir[2],-dir[2],dir[1])
end

@recipe function f(b::Ellipse)
    grid := false
    axis := false
    ticks := false
    legend := false
    aspect_ratio := :equal
    linewidth --> 0.2
    linecolor --> :black
    fillcolor --> :black
    fillalpha --> 0.5
    fillrange := 0.0
    @series begin
        t = 0:0.1:2π
        bds = b.rotate*[b.radii[1]*cos.(t)'; b.radii[2]*sin.(t)'] .+ b.center
        bds = hcat(bds, bds[:, 1])
        bds[1,:], bds[2,:]
    end
end

end
