module SymbolicCurves
using CurveProximityQueries
using IntervalArithmetic
using ForceImport
using StaticArrays
@force using Reduce.Algebra

export SymCurve

struct SymCurve{D, T} <: Curve{D, T}
    f
    s
    limits::Interval{T}
end

function SymCurve(xexpr, yexpr, limits::Interval{T}, θ, xcenter, ycenter) where {T}
    dx = df(xexpr, :t)
    dy = df(yexpr, :t)
    sexpr = int(dx*dx + dy*dy, :t)
    s = eval(:(t->$sexpr))
    x = eval(:(t->$xexpr))
    y = eval(:(t->$yexpr))
    rx(u) =  cos(θ)*x(u) + sin(θ)*y(u) + xcenter
    ry(u) = -sin(θ)*x(u) + cos(θ)*y(u) + ycenter
    F = SymCurveBase(rx, ry)
    S = SymArclengthCurve(s)
    SymCurve{2, T}(F, S, limits)
end


struct SymCurveBase
    x
    y
end

struct SymArclengthCurve
    g
end

function (b::SymCurveBase)(t::Real)
    SVector{2, Float64}(b.x(t), b.y(t))
end

(b::SymArclengthCurve)(t) = b.g(t)

(b::SymCurve)(t) = b.f(t)

import CurveProximityQueries: arclength
function arclength(b::SymCurve, θ::Interval)
    s = b.s(θ.hi) - b.s(θ.lo)
    s > 0 ? sqrt(diam(θ)*s) : 0.0
end

#  function ranunculoid(a::Real, θ, xcenter, ycenter)
#      @vars t
#      x(t) = a*0.05*(6*cos(t) - cos(6*t))
#      y(t) = a*0.05*(6*sin(t) - sin(6*t))
#      SymCurve(x, y, Interval(0, 2π), θ, xcenter, ycenter)
#  end
#
#  function fishcurve(a::Real, θ, xcenter, ycenter)
#      @vars t
#      x(t) = a*0.3*(cos(t) - sin(t)^2/√2)
#      y(t) = a*0.3*sin(t)*cos(t)
#      SymCurve(x, y, Interval(0, 2π), θ, xcenter, ycenter)
#  end
#
#  function triplesine(a::Real, b::Real, c::Real, θ, xcenter, ycenter)
#      @vars t
#      x(t) = 0.2*(cos(a*t) + cos(b*t)/2 + sin(c*t)/3)
#      y(t) = 0.2*(sin(a*t) + sin(b*t)/2 + cos(c*t)/3)
#      SymCurve(x, y, Interval(0, 2π), θ, xcenter, ycenter)
#  end
#
#  function lissajous(p::Real, q::Real, θ, xcenter, ycenter)
#      @vars t
#      x(t) = 0.3*sin(p*t)
#      y(t) = 0.3*sin(q*t)
#      SymCurve(x, y, Interval(0, 2π), θ, xcenter, ycenter)
#  end
#
#  function clothoid(a::Real, θ, xcenter, ycenter)
#      x(u) = a*(sqrt(2) * sqrt(pi) * fresnelc((sqrt(2) * u) / sqrt(pi)) * gamma(1 / 4)) / (8 * gamma(5 / 4))
#      y(u) = a*(3 * sqrt(2) * sqrt(pi) * fresnels((sqrt(2) * u) / sqrt(pi)) * gamma(3 / 4)) / (8 * gamma(7 / 4))
#      s(u) = a*a*u
#      rx(u) =  cos(θ)*x(u) + sin(θ)*y(u) + xcenter
#      ry(u) = -sin(θ)*x(u) + cos(θ)*y(u) + ycenter
#      F = SymCurveBase(rx, ry)
#      S = SymArclengthCurve(s)
#      SymCurve(F, S, Interval(-2π, 2π))
#  end
end
