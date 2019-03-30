import Base: show, product, @_inline_meta
import Base: eltype, getindex
import Base: +,-,*,/,sum

struct Bernstein{D, N, T} <: ParametricCurve{D, T}
    control_points::SVector{N, SVector{D, T}}
end

function Bernstein(cp::Vararg{SVector,N}) where {N}
    cpts = SVector{N}(cp...)
    Bernstein(cpts)
end

Base.eltype(::Type{Bernstein{D, N, T}}) where {D, N, T} = T
Base.getindex(b::Bernstein, j::Int) = b.control_points[j]

function Base.show(io::IO, b::Bernstein{D, N, T}) where {D, N, T}
    ordind = (N-1)%10 == 1 ? "st" : "th"
    ordind = (N-2)%10 == 2 ? "nd" : ordind
    ordind = (N-3)%10 == 3 ? "rd" : ordind
    ordind = 11 ≤ (N-1)%100 ≤ 13 ? "th" : ordind

    print(io, "a ", (N-1), ordind, " order Bernstein polynomial with control points at:\n",
          b.control_points.data)
end

function _eval(b::Bernstein{D, N}, t::Real) where {D, N}
    h = 1.0
    u = 1.0 - t
    q = b[1]
    if t < 0.5
        u = t/u
        for k = 1:N-1
            h = h*u*(N-k)
            h = h/(k + h)
            q = (1.0-h)*q + h*b[k+1]
        end
    else
        u = u/t
        for k = 1:N-1
            h = h*(N-k)
            h = h/(k*u + h)
            q = (1.0-h)*q + h*b[k+1]
        end
    end
    return q
end

function (b::Bernstein)(t::Real)
    _eval(b, t)
end

function (b::Bernstein{1})(t::Real)
    _eval(b, t)[1]
end

@generated function differentiate(b::Bernstein{D, N, T}) where {D, N, T}
    exprs = Array{Expr}(undef, N-1)
    for k = 1:N-1
        exprs[k] = :((b[$(k+1)] - b[$k])*(N-1))
    end

    return quote
        @_inline_meta
        @inbounds elements = tuple($(exprs...))
        @inbounds return Bernstein(SVector{N-1}(elements))
    end
end

@generated function integrate(b::Bernstein{D, N, T}) where {D, N, T}
    exprs = Array{Expr}(undef, N+1)
    for k = 1:N+1
        exprs[k] = :(zero(b[1]))
        for j=1:k-1
            exprs[k] = :($(exprs[k]) + b[$j])
        end
        exprs[k] = :($(exprs[k])/N)
    end

    return quote
        @_inline_meta
        @inbounds elements = tuple($(exprs...))
        @inbounds return Bernstein(SVector{N+1}(elements))
    end
end

function *(k::Real, b::Bernstein)
    Bernstein(b.control_points*k)
end
*(b::Bernstein, k::Real) = k*b
/(b::Bernstein, k::Real) = (1/k)*b
-(b::Bernstein) = (-1)*b
+(b::Bernstein) = (+1)*b

function +(a::Bernstein{D, Na}, b::Bernstein{D, Nb}) where {D, Na, Nb}
    if Na < Nb
        a = _degelevate(a)
    elseif Na > Nb
        b = _degelevate(b)
    end
    a + b
end

function +(a::Bernstein{D, N}, b::Bernstein{D, N}) where {D, N}
    Bernstein(a.control_points + b.control_points)
end

function -(a::Bernstein{D}, b::Bernstein{D}) where {D}
    a + -b
end

@generated function _degelevate(b::Bernstein{D, N, T}) where {D, N, T}
    exprs = Array{Expr}(undef, N+1)
    exprs[1] = :(b[1])
    exprs[N+1] = :(b[N])
    for k = 2:N
        exprs[k] = :($(k/N)*b[$(k-1)] + (1 - $(k/N))*b[$k])
    end

    return quote
        @_inline_meta
        @inbounds elements = tuple($(exprs...))
        @inbounds return Bernstein(SVector{N+1}(elements))
    end
end

@generated function _binomial(::Type{Bernstein{D, N, T}}) where {D, N, T}
    nck = Array{Expr}(undef, N)
    for i = 1:N
        nck[i] = (i==1) ? :(one($T)) : :($(nck[i-1])*($N-$i+1)/($i-1))
    end
    return nck
end

@generated function _magsquared(b::Bernstein{D, N, T}) where {D, N, T}
    exprs = Array{Expr}(undef, 2N-1)
    mck = Array{Expr}(undef, N)
    den = :nothing

    nck = _binomial(b)

    for k = 1:2N-1
        exprs[k] = :(zero($T))
        den = (k==1) ? :(one($T)) : :($den*$(2N-k)/$(k-1))
        for j = max(1, k-N+1):min(N, k)
            exprs[k] = :($(exprs[k]) + $(nck[j])*$(nck[k-j+1])*b[$j]'*b[$(k-j+1)])
        end
        exprs[k] = :(SVector{1}($(exprs[k])/$den))
    end

    return quote
        @_inline_meta
        @inbounds elements = tuple($(exprs...))
        @inbounds return Bernstein(SVector{2N-1}(elements))
    end
end

struct RectifiableBernstein{D, N, M, T} <: RectifiableCurve{D, T}
    f::Bernstein{D, N, T}
    s::Bernstein{1, M, T}
    limits::Interval{T}
end

function RectifiableBernstein(cp::Vararg{SVector,N}; limits=Interval(0., 1.)) where {N}
    cpts = SVector{N}(cp...)
    f = Bernstein(cpts)
    df = differentiate(f)
    s = integrate(_magsquared(df))
    RectifiableBernstein(f, s, limits)
end

function RectifiableBernstein(f::Bernstein; limits=Interval(0., 1.)) where {N}
    df = differentiate(f)
    s = integrate(_magsquared(df))
    RectifiableBernstein(f, s, limits)
end

scalet(b::RectifiableBernstein, t) = (t - b.limits.lo)/diam(b.limits)

function (b::RectifiableBernstein)(t::Real)
    t = scalet(b, t)
    b.f(t)
end

function arclength(b::RectifiableBernstein, t::Real)
    t = scalet(b, t)
    s = b.s(t)
    s > 0 ? sqrt(t*s) : 0.0
end

function arclength(b::RectifiableBernstein, θ::Interval)
    θ = scalet(b, θ)
    s = b.s(θ.hi) - b.s(θ.lo)
    s > 0 ? sqrt(diam(θ)*s) : 0.0
end


Base.eltype(::Type{RectifiableBernstein{D, N, T}}) where {D, N, T} = T
Base.getindex(b::RectifiableBernstein, j::Int) = b.f.control_points[j]

function Base.show(io::IO, b::RectifiableBernstein{D, N, T}) where {D, N, T}
    ordind = (N-1)%10 == 1 ? "st" : "th"
    ordind = (N-2)%10 == 2 ? "nd" : ordind
    ordind = (N-3)%10 == 3 ? "rd" : ordind
    ordind = 11 ≤ (N-1)%100 ≤ 13 ? "th" : ordind

    print(io, "a ", (N-1), ordind, " order Rectifiable Bernstein polynomial with control points at:\n",
          b.f.control_points.data, "\nwith an arclength of ", arclength(b, b.limits.hi))
end

function *(k::Real, b::RectifiableBernstein)
    RectifiableBernstein(b.f*k, limits=b.limits)
end
*(b::RectifiableBernstein, k::Real) = k*b
/(b::RectifiableBernstein, k::Real) = (1/k)*b
-(b::RectifiableBernstein) = (-1)*b
+(b::RectifiableBernstein) = (+1)*b

function +(a::RectifiableBernstein{D}, b::RectifiableBernstein{D}) where {D}
    @assert a.limits == b.limits
    c = a.f + b.f
    RectifiableBernstein(c, limits=a.limits)
end

function -(a::RectifiableBernstein{D}, b::RectifiableBernstein{D}) where {D}
    a + -b
end

import Random: rand
function rand(::Type{Bernstein{D, N}}) where {D, N}
    pts = [@SVector(rand(D)) for i=1:N]
    Bernstein(pts...)
end

function rand(::Type{RectifiableBernstein{D, N}}) where {D, N}
    bern = rand(Bernstein{D, N})
    RectifiableBernstein(bern)
end
