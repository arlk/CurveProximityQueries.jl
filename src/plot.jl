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
        xy.(θ)
    end
end

@recipe function f(a::ConvexPolygon)
    grid --> false
    axis --> false
    ticks --> false
    legend --> false
    aspect_ratio --> :equal
    linewidth --> 1.0
    linecolor --> :black
    fillcolor --> :black
    fillalpha --> 0.5
    fillrange := 0.0
    leg := false
    @series begin
        vertices = vcat(a.pts, a.pts[1:1])
        vertices
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
    markerstrokewidth --> 0.0
    @series begin
        v = [a...]
        v
    end
end
