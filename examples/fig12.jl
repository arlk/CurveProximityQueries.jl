using BenchmarkTools
using ConvexBodyProximityQueries
using CurveProximityQueries
using Plots
using StaticArrays
import Random: seed!

pgfplots()

function randomtraj(num)
    Bs = []
    xax = collect(0.0:0.2:1.0)'
    i1 = @SVector [-0.2, 0.5]
    i2 = @SVector [ 0.0, 0.5]
    e1 = @SVector [ 1.0, 0.5]
    e2 = @SVector [ 1.2, 0.5]
    for i = 1:num
        m1 = @SVector rand(2)
        m2 = @SVector rand(2)
        b = Bernstein(i1, i2, m1, m2, e1, e2)
        push!(Bs, b)
    end
    return Bs
end

import CurveProximityQueries:tolerance_verification
function tolerance_verification(obs::Array, traj::Curve, Δ::Real)
    ret = true
    for o in obs
        ret = ret && tolerance_verification(o, traj, Δ)
    end
    return ret
end

import CurveProximityQueries:collision_detection
function collision_detection(obs::Array, traj::Curve)
    ret = false
    for o in obs
        ret = ret || collision_detection(o, traj)
    end
    return ret
end

function counttv(obs, trajs)
    cnt = 0.0
    for t in trajs
        cnt += tolerance_verification(obs, t, 0.03)
    end
    return cnt
end

seed!(22)
obs = [randpoly([0.35, 0.6]; n=20, scale=0.2),
       randpoly([0.65, 0.4]; n=20, scale=0.2)]
trajs = randomtraj(1000)

#  @btime counttv($obs, $trajs)

plot()
plot!.(obs,
       fillalpha=1.0,
       fillcolor=Gray(0.7),
      )
cnt = zeros(3)
for t in trajs
    if tolerance_verification(obs, t, 0.03)
        cnt[1] += 1
        plot!(t,
              linecolor = RGBA(77/255,141/225,205/255,1.0),
              linewidth = 1.2)
    elseif !collision_detection(obs, t)
        cnt[2] += 1
        plot!(t,
              linecolor = RGBA(231/255,95/225,97/255,0.3),
              linewidth = 1.2)
    else
        cnt[3] += 1
        plot!(t,
              linecolor = GrayA(0.5, 0.05),
              linewidth = 0.5)
    end
end
plot!(xlim=(-0.3, 1.3), ylim=(0.0, 1.0))
savefig("out/fig12.pdf")
