using RecipesBase

@recipe function f(b::Curve)
    grid --> false
    axis --> false
    ticks --> false
    legend --> false
    aspect_ratio --> :equal
    linewidth --> 2.0
    @series begin
        xy(u) = b(u)
        lim = b.limits
        Δ = diam(lim)/100
        θ = lim.lo:Δ:lim.hi
        Tuple.(xy.(θ))
    end
end

@recipe function f(b::Curve, lim::Interval)
    grid --> false
    axis --> false
    ticks --> false
    legend --> false
    aspect_ratio --> :equal
    linewidth --> 2.0
    @series begin
        xy(u) = b(u)
        Δ = diam(lim)/100
        θ = lim.lo:Δ:lim.hi
        Tuple.(xy.(θ))
    end
end

@recipe function f(a::Tuple{SVector{D, T}, SVector{D, T}}) where {D, T}
    grid --> false
    axis --> false
    ticks --> false
    legend --> false
    aspect_ratio --> :equal
    linewidth --> 1.0
    linecolor --> :red
    markercolor --> :red
    markershape --> :circle
    markersize --> 5.0
    markerstrokecolor --> :red
    @series begin
        Tuple.([a...])
    end
end
