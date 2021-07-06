using CurveProximityQueries
using ConvexBodyProximityQueries
using Test
using LinearAlgebra
using StaticArrays
using IntervalArithmetic
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
        obs = @point zeros(2)
        c = Bernstein([ [0.23603334566204692, 0.34651701419196046],
                        [0.3127069683360675, 0.00790928339056074],
                        [0.4886128300795012, 0.21096820215853596],
                        [0.951916339835734, 0.9999046588986136],
                        [0.25166218303197185, 0.9866663668987996],
                        [0.5557510873245723, 0.43710797460962514],
                        [0.42471785049513144, 0.773223048457377] ])
        @test abs(minimum_distance(obs, c) - 0.3683522584741768) ≤ 1e-5
        @test tolerance_verification(obs, c, 0.3) == true
        @test tolerance_verification(obs, c, 0.4) == false
        @test collision_detection(obs, c) == false
        @test abs(minimum_distance(c, obs) - 0.3683522584741768) ≤ 1e-5
        @test tolerance_verification(c, obs, 0.3) == true
        @test tolerance_verification(c, obs, 0.4) == false
        @test collision_detection(c, obs) == false

        obs = @point zeros(3)
        c = Bernstein([ [0.23603334566204692, 0.34651701419196046, 0.3127069683360675],
                        [0.00790928339056074, 0.4886128300795012, 0.21096820215853596],
                        [0.951916339835734, 0.9999046588986136, 0.25166218303197185],
                        [0.9866663668987996, 0.5557510873245723, 0.43710797460962514],
                        [0.42471785049513144, 0.773223048457377, 0.2811902322857298],
                        [0.20947237319807077, 0.25137920979222494, 0.02037486871266725],
                        [0.2877015122756894, 0.859512136087661, 0.07695088688120899],
                        [0.6403962459899388, 0.8735441302706854, 0.27858242002877853],
                        [0.7513126327861701, 0.6448833539420931, 0.07782644396003469],
                        [0.8481854810000327, 0.0856351682044918, 0.5532055454580578],
                        [0.46335024592359875, 0.18582130997265378, 0.11198087695816716] ])
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
        c = Bernstein([ [0.23603334566204692, 0.34651701419196046],
                        [0.3127069683360675, 0.00790928339056074],
                        [0.4886128300795012, 0.21096820215853596],
                        [0.951916339835734, 0.9999046588986136],
                        [0.25166218303197185, 0.9866663668987996],
                        [0.5557510873245723, 0.43710797460962514],
                        [0.42471785049513144, 0.773223048457377] ])
        d = Bernstein([ [0.2811902322857298, 0.20947237319807077],
                        [0.25137920979222494, 0.02037486871266725],
                        [0.2877015122756894, 0.859512136087661] ])
        @test abs(minimum_distance(c, d)) ≤ 1e-5
        @test tolerance_verification(c, d, 0.3) == false
        @test collision_detection(c, d) == true

        obs = @point zeros(3)
        c = Bernstein([ [0.23603334566204692, 0.34651701419196046, 0.3127069683360675],
                        [0.00790928339056074, 0.4886128300795012, 0.21096820215853596],
                        [0.951916339835734, 0.9999046588986136, 0.25166218303197185],
                        [0.9866663668987996, 0.5557510873245723, 0.43710797460962514],
                        [0.42471785049513144, 0.773223048457377, 0.2811902322857298],
                        [0.20947237319807077, 0.25137920979222494, 0.02037486871266725],
                        [0.2877015122756894, 0.859512136087661, 0.07695088688120899],
                        [0.6403962459899388, 0.8735441302706854, 0.27858242002877853],
                        [0.7513126327861701, 0.6448833539420931, 0.07782644396003469],
                        [0.8481854810000327, 0.0856351682044918, 0.5532055454580578],
                        [0.46335024592359875, 0.18582130997265378, 0.11198087695816716] ])
        d = Bernstein([ [0.976311881619359, 0.051614620674327094, 0.5380295812064833],
                        [0.4556920516275036, 0.2793951106725605, 0.17824610354168602],
                        [0.5489828719625274, 0.37097066286146885, 0.8941659192657594],
                        [0.6480537482231894, 0.4170393538841062, 0.14456554241360564],
                        [0.6224031828206811, 0.8723344353741975, 0.5249746566167794] ])
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
