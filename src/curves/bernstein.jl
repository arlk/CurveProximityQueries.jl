import Base: show, product, @_inline_meta
import Base: eltype, getindex
import Base: +,-,*,/,sum

struct BernsteinBase{D, N, T, A <: StaticVector{D, T}}
    control_points::SVector{N, A}
end

function BernsteinBase(cp::Vararg{StaticVector{<:Any, <:Number}})
    cpts = SVector(map(t->map(float, t), cp)...)
    BernsteinBase(cpts)
end

function BernsteinBase(data::AbstractVector)
    @assert !isempty(data)
    cptsz = length(data)
    dimsz = length(first(data))
    T1 = eltype(data)
    T2 = float(eltype(T1))
    T = StaticArrays.similar_type(T1, T2, StaticArrays.Size((dimsz,)))
    cpts = SVector{cptsz}((T(d) for d in data)...)
    BernsteinBase(cpts)
end

Base.eltype(::Type{BernsteinBase{D, N, T}}) where {D, N, T} = T
Base.getindex(b::BernsteinBase, j::Int) = b.control_points[j]

function _eval(b::BernsteinBase{D, N}, t::Real) where {D, N}
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

function (b::BernsteinBase)(t::Real)
    _eval(b, t)
end

function (b::BernsteinBase{1})(t::Real)
    _eval(b, t)[1]
end

@generated function differentiate(b::BernsteinBase{D, N, T}) where {D, N, T}
    exprs = Array{Expr}(undef, N-1)
    for k = 1:N-1
        exprs[k] = :((b[$(k+1)] - b[$k])*(N-1))
    end

    return quote
        @_inline_meta
        @inbounds elements = tuple($(exprs...))
        @inbounds return BernsteinBase(SVector{N-1}(elements))
    end
end

@generated function integrate(b::BernsteinBase{D, N, T}) where {D, N, T}
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
        @inbounds return BernsteinBase(SVector{N+1}(elements))
    end
end

function *(k::Real, b::BernsteinBase)
    BernsteinBase(b.control_points*k)
end
*(b::BernsteinBase, k::Real) = k*b
/(b::BernsteinBase, k::Real) = (1/k)*b
-(b::BernsteinBase) = (-1)*b
+(b::BernsteinBase) = (+1)*b

function +(a::BernsteinBase{D, Na}, b::BernsteinBase{D, Nb}) where {D, Na, Nb}
    if Na < Nb
        a = _degelevate(a)
    elseif Na > Nb
        b = _degelevate(b)
    end
    a + b
end

function +(a::BernsteinBase{D, N}, b::BernsteinBase{D, N}) where {D, N}
    BernsteinBase(a.control_points + b.control_points)
end

function -(a::BernsteinBase{D}, b::BernsteinBase{D}) where {D}
    a + -b
end

@generated function _degelevate(b::BernsteinBase{D, N, T}) where {D, N, T}
    exprs = Array{Expr}(undef, N+1)
    exprs[1] = :(b[1])
    exprs[N+1] = :(b[N])
    for k = 2:N
        exprs[k] = :($(k/N)*b[$(k-1)] + (1 - $(k/N))*b[$k])
    end

    return quote
        @_inline_meta
        @inbounds elements = tuple($(exprs...))
        @inbounds return BernsteinBase(SVector{N+1}(elements))
    end
end

@generated function _binomial(::Type{BernsteinBase{D, N, T, A}}) where {D, N, T, A}
    nck = Array{Expr}(undef, N)
    for i = 1:N
        nck[i] = (i==1) ? :(one($T)) : :($(nck[i-1])*($N-$i+1)/($i-1))
    end
    return nck
end

@generated function _magsquared(b::BernsteinBase{D, N, T}) where {D, N, T}
    exprs = Array{Expr}(undef, 2N-1)
    mck = Array{Expr}(undef, N)
    den = :nothing

    nck = _binomial(b)

    for k = 1:2N-1
        exprs[k] = :(zero($T)*zero($T))
        den = (k==1) ? :(one($T)) : :($den*$(2N-k)/$(k-1))
        for j = max(1, k-N+1):min(N, k)
            exprs[k] = :($(exprs[k]) + $(nck[j])*$(nck[k-j+1])*b[$j]'*b[$(k-j+1)])
        end
        exprs[k] = :(SVector{1}($(exprs[k])/$den))
    end

    return quote
        @_inline_meta
        @inbounds elements = tuple($(exprs...))
        @inbounds return BernsteinBase(SVector{2N-1}(elements))
    end
end

struct Bernstein{D, N, M, Tf, Ts, Af, As, J} <: Curve{D, Tf}
    f::BernsteinBase{D, N, Tf, Af}
    s::BernsteinBase{1, M, Ts, As}
    limits::Interval{J}
end

# <:Number specificity needed to avoid error when doing
# Bernstein(@SVector([ @SVector([a,b]), @SVector([c,d]), ... ]))
function Bernstein(cp::Vararg{StaticVector{<:Any, <:Number}}; limits=Interval(0., 1.))
    Bernstein(BernsteinBase(cp...), limits=limits)
end

function Bernstein(f::BernsteinBase; limits=Interval(0., 1.))
    df = differentiate(f)
    s = integrate(_magsquared(df))
    Bernstein(f, s, limits)
end

function Bernstein(data::AbstractVector; limits=Interval(0., 1.))
    Bernstein(BernsteinBase(data), limits=limits)
end

scalet(b::Bernstein, t) = (t - b.limits.lo)/diam(b.limits)

function (b::Bernstein)(t::Real)
    t = scalet(b, t)
    b.f(t)
end

function arclength(b::Bernstein, t::Real)
    t = scalet(b, t)
    s = b.s(t)
    s > zero(s) ? sqrt(t*s) : float(zero(eltype(b)))
end

function arclength(b::Bernstein, θ::Interval)
    θ = scalet(b, θ)
    s = b.s(θ.hi) - b.s(θ.lo)
    s > zero(s) ? sqrt(diam(θ)*s) : float(zero(eltype(b)))
end

function arclength(b::Bernstein)
    s = b.s(1.0)
    s > zero(s) ? sqrt(s) : float(zero(eltype(b)))
end

Base.eltype(::Type{<:Bernstein{D, N, M, T}}) where {D, N, M, T} = T
Base.getindex(b::Bernstein, j::Int) = b.f.control_points[j]

function Base.show(io::IO, b::Bernstein{D, N}) where {D, N}
    ordind = (N-1)%10 == 1 ? "st" : "th"
    ordind = (N-1)%10 == 2 ? "nd" : ordind
    ordind = (N-1)%10 == 3 ? "rd" : ordind
    ordind = 11 ≤ (N-1)%100 ≤ 13 ? "th" : ordind

    print(io, "a ", (N-1), ordind, " order Bernstein polynomial with control points at:\n",
          b.f.control_points.data, "\nwith an arclength of ", arclength(b, b.limits.hi))
end

function *(k::Real, b::Bernstein)
    Bernstein(b.f*k, limits=b.limits)
end
*(b::Bernstein, k::Real) = k*b
/(b::Bernstein, k::Real) = (1/k)*b
-(b::Bernstein) = (-1)*b
+(b::Bernstein) = (+1)*b

function +(a::Bernstein{D}, b::Bernstein{D}) where {D}
    @assert a.limits == b.limits
    c = a.f + b.f
    Bernstein(c, limits=a.limits)
end

function -(a::Bernstein{D}, b::Bernstein{D}) where {D}
    a + -b
end

function differentiate(b::Bernstein)
    c = differentiate(b.f)/diam(b.limits)
    Bernstein(c, limits=b.limits)
end

import Random: rand
function rand(::Type{BernsteinBase{D, N}}) where {D, N}
    pts = [@SVector(rand(D)) for i=1:N]
    BernsteinBase(pts...)
end

function rand(::Type{Bernstein{D, N}}; limits=Interval(0., 1.)) where {D, N}
    bern = rand(BernsteinBase{D, N})
    Bernstein(bern, limits=limits)
end
