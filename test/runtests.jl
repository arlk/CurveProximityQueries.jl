using CurveProximityQueries
using ConvexBodyProximityQueries
using Test
using LinearAlgebra
using StaticArrays
using IntervalArithmetic
using Random: seed!
using Unitful: nm, μm, mm, cm

import CurveProximityQueries: differentiate, integrate

@testset "CurveProximityQueries" begin
    @testset "Constructors" begin
        # should return SVector when constructed without SVectors
        @test Bernstein([[0.0, 0.0], [1.0, 1.0]])(0.5) === @SVector([0.5, 0.5])

        # non-floats arguments should become floats to satisfy assumptions
        # elsewhere in the code (e.g. use of `eps(T)`)
        @test eltype(Bernstein([[0, 0], [1, 1]])) <: Float64

        # none of the following should error; some should be inferrable
        @test Bernstein([
                @SVector([0.0, 0.0]),
                @SVector([1.0, 1.0])
            ])(0.5) === @SVector([0.5, 0.5])
        @test @inferred(Bernstein(@SVector([
                @SVector([0.0, 0.0]),
                @SVector([1.0, 1.0])
            ])))(0.5) === @SVector([0.5, 0.5])
        @test @inferred(Bernstein(
                @SVector([0.0, 0.0]),
                @SVector([1.0, 1.0])
            ))(0.5) === @SVector([0.5, 0.5])

        # similar reasoning below, but allow for any StaticVector
        struct Point{T} <: StaticArrays.FieldVector{2, T}
            x::T
            y::T
        end
        StaticArrays.similar_type(::Type{P}, ::Type{T},
            ::StaticArrays.Size{(2,)}) where {P <: Point, T} = Point{T}

        @test Bernstein([Point(0.0, 0.0), Point(1.0, 1.0)])(0.5) ===
            Point(0.5, 0.5)
        @test @inferred(
            Bernstein(@SVector([Point(0.0, 0.0), Point(1.0, 1.0)])))(0.5) ===
                Point(0.5, 0.5)
        @test @inferred(Bernstein(Point(0.0, 0.0), Point(1.0, 1.0)))(0.5) ===
            Point(0.5, 0.5)
    end
    @testset "Bernstein Polynomials" begin
        B2 = rand(Bernstein{2,8})
        B3 = rand(Bernstein{3,5})
        @test typeof(B2) <: Bernstein{2,8}
        @test typeof(B3) <: Bernstein{3,5}
        @test sum(norm.((differentiate(integrate(B2.f)) - B2.f).control_points)) ≤ 1e-10
        @test sum(norm.((differentiate(integrate(B3.f)) - B3.f).control_points)) ≤ 1e-10
    end
    @testset "Curve - Polygon" begin
        seed!(1)
        obs = @point zeros(2)
        c = rand(Bernstein{2,7})
        @test abs(minimum_distance(obs, c) - 0.3683522584741768) ≤ 1e-5
        @test tolerance_verification(obs, c, 0.3) == true
        @test tolerance_verification(obs, c, 0.4) == false
        @test collision_detection(obs, c) == false
        @test abs(minimum_distance(c, obs) - 0.3683522584741768) ≤ 1e-5
        @test tolerance_verification(c, obs, 0.3) == true
        @test tolerance_verification(c, obs, 0.4) == false
        @test collision_detection(c, obs) == false

        seed!(1)
        obs = @point zeros(3)
        c = rand(Bernstein{3,11})
        @test abs(minimum_distance(obs, c) - 0.511627527056288) ≤ 1e-5
        @test tolerance_verification(obs, c, 0.5) == true
        @test tolerance_verification(obs, c, 0.6) == false
        @test collision_detection(obs, c) == false
        @test abs(minimum_distance(c, obs) - 0.511627527056288) ≤ 1e-5
        @test tolerance_verification(c, obs, 0.5) == true
        @test tolerance_verification(c, obs, 0.6) == false
        @test collision_detection(c, obs) == false
    end
    @testset "Curve - Curve" begin
        seed!(1)
        c = rand(Bernstein{2,7})
        d = rand(Bernstein{2,3})
        @test abs(minimum_distance(c, d)) ≤ 1e-5
        @test tolerance_verification(c, d, 0.3) == false
        @test collision_detection(c, d) == true

        seed!(1)
        obs = @point zeros(3)
        c = rand(Bernstein{3,11})
        d = rand(Bernstein{3,5})
        @test abs(minimum_distance(c, d) - 0.137873546917387) ≤ 1e-5
        @test tolerance_verification(c, d, 0.1) == true
        @test tolerance_verification(c, d, 0.3) == false
        @test collision_detection(c, d) == false
    end
    @testset "Dimensionful curves" begin
        b2 = Bernstein([[0.0, 1.0]μm, [1.0, 1.0]μm, [1.0, 0.0]μm])
        @test b2(0.0) == [0, 1000]nm
        @test b2(1.0) == [0.001, 0]mm
        @test arclength(b2) == arclength(Bernstein([[0,1], [1,1], [1,0]])) * μm
        @test length(CurveProximityQueries.cvxhull(b2, Interval(0.0, 1.0))) == 8

        b3 = Bernstein([[0.0, 0.0, 0.0]mm, [1.0, 0.0, 0.0]mm, [1.0, 0.0, 1.0]mm])
        c3 = Bernstein([[2.0, 1.0, 1.0]cm, [2.0, 3.0, 2.0]cm, [3.0, 1.0, 1.0]cm])
        @test length(CurveProximityQueries.cvxhull(b3, Interval(0.0, 1.0))) == 16
        @test abs(minimum_distance(b3, c3, atol=1e-8mm) - 2.328cm) ≤ 1e-4cm
        @test tolerance_verification(b3, c3, 2cm, atol=1e-8mm) == true
        @test collision_detection(b3, c3, atol=1e-8mm) == false
    end
    @testset "Support" begin
        @test CurveProximityQueries.support(Point(1.0, 2.0), Point(1.0, 1.0)) ===
            Point(1.0, 2.0)
        @test CurveProximityQueries.support(@SVector(
            [Point(0.0, 1.0), Point(1.0, 2.0)]), @SVector([1.0, 1.0])) ===
                Point(1.0, 2.0)
        @test CurveProximityQueries.support(@SVector(
            [@SVector([0.0, 1.0]), @SVector([1.0, 2.0])]),
                @SVector([1.0, 1.0])) === @SVector([1.0, 2.0])
        @test CurveProximityQueries.support(@SVector(
            [Point(0.0, 1.0), Point(1.0, 2.0), Point(3.0, 4.0)]), Point(1.0, 1.0)) ===
                Point(3.0, 4.0)

    end
end
