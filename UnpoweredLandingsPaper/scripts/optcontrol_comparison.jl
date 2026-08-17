## compare different methods to simulate optimal solutions
using UnpoweredLandingsPaper, PointMass, BoundaryValueDiffEq, OrdinaryDiffEqTsit5, NonlinearSolve, CairoMakie
import CairoMakie.Makie.SpecApi as S
using ForwardDiff: derivative, gradient

struct OptControlAug <: Control end

function PointMass.alpha(x, p::Parameters{<:Aerodynamics, <:OptControlAug}, t)
    ν, γ, _, _, p₁, p₂ = x
    f(α, p) = derivative(α -> (1 + p₁) * (-CD(α, p) * ν^2 - sin(γ)) + p₂ * (CL(α, p) * ν - cos(γ) / ν), α)
    prob = IntervalNonlinearProblem(f, eltype(x)[0.001, π / 2 + 0.001], p)
    return solve(prob, ITP()).u
end
PointMass.thrust(x, p::Parameters{<:Aerodynamics, <:OptControlAug}, t) = 0.0
PointMass.thrust_angle(x, p::Parameters{<:Aerodynamics, <:OptControlAug}, t) = 0.0

H(x, p, t) = (1 + x[5]) * ν̇(x, p, t) + x[6] * γ̇(x, p, t)

function full_dynamics!(dx, x, p::Parameters{<:Aerodynamics, OptControlAug}, t)
    ∇x = gradient(x -> H(x, p, t), x)
    dx[1] = ∇x[5]
    dx[2] = ∇x[6]
    dx[3] = x[1] * cos(x[2])
    dx[4] = x[1] * sin(x[2])
    dx[5] = -∇x[1]
    dx[6] = -∇x[2]
    return nothing
end


function optcontrol_comparison()
    ## parameters and boundary conditions
    pa = TrigAerodynamics(1.5, 0.03, 2.0)
    p1 = Parameters(pa, OptControlAug())
    p2 = Parameters(pa, OptControl())
    cb = ContinuousCallback((u, t, int) -> u[2], nothing, int -> terminate!(int))
    tspan = (0.0, -3.0)
    νtd = 0.1
    γtd = asin(-max_CD(pa) * νtd^2)
    u₁ = [νtd, γtd, 0.0, 0.0, 0.0, 0.0] # ν, γ, χ, ζ, p₁, p₂

    ## solve augmented system backwards as IVP from touchdown
    sol_aug = simulate(full_dynamics!, u₁, p1, tspan; callback = cb)

    ## solve original system backwards as IVP from touchdown solving
    sol_org = simulate(PointMass.full_dynamics!, u₁, p2, tspan; callback = cb)

    ## solve with direct collocation
    sol_drc = solve_ocp(dc_min_speed_unconstrained(sol_org.u[1], pa); solver = :madnlp, grid_size = 200, display = false)

    # plot
    f = Figure()
    sols = [sol_aug, sol_org, sol_drc]
    labels = ["BVP", "IVP", "DC"]
    ls = [:solid, :dash, :dot]
    colors = cs1
    ax1 = trajectory_axis(
        mapreduce(vcat, sols, colors, ls, labels) do sol, color, linestyle, label
            trajectory(sol, 0; color, linestyle, label)
        end
    )

    ax2a = timeseries_axis(
        map(sols, colors, ls, labels) do sol, color, linestyle, label
            pts = [Point2(t, u[1]) for (u, t) in zip(sol.u, sol.t)]
            S.Lines(pts; color, linestyle, label)
        end, "ν"
    )
    ax2b = timeseries_axis(
        map(sols, colors, ls, labels) do sol, color, linestyle, label
            pts = [Point2(t, u[2]) for (u, t) in zip(sol.u, sol.t)]
            S.Lines(pts; color, linestyle, label)
        end, "γ"
    )
    ax2c = timeseries_axis(
        map(sols, colors, ls, labels) do sol, color, linestyle, label
            _f = sol isa ODESolution ? t -> PointMass.alpha(sol(t), sol.prob.p, t) : sol.cntrlf
            pts = [Point2f(t, _f(t)) for (u, t) in zip(sol.u, sol.t)]
            S.Lines(pts; color, linestyle, label)
        end, "α"
    )

    plot(
        f[1, 1], S.GridLayout(
            (1, 1) => ax1,
            (1, 2) => S.GridLayout([ax2a, ax2b, ax2c]),
        )
    )

    Legend(f[2, 1], f.content[1]; orientation = :horizontal)

    return f
end


# ## solve augmented system as BVP
# tspan = extrema(sol_aug.t)
# uguess = sol_aug.u[1]

# function bca!(resid_a, u_a, p)
#     resid_a[1] = u_a[1] - sol_aug.u[1][1]
#     resid_a[2] = u_a[2] - sol_aug.u[1][2]
#     return nothing
# end
# function bcb!(resid_b, u_b, p)
#     resid_b[1] = u_b[3]
#     resid_b[2] = u_b[4]
#     resid_b[3] = u_b[5]
#     resid_b[4] = u_b[6]
#     return nothing
# end
# prob = TwoPointBVProblem(
#     full_dynamics!, (bca!, bcb!), uguess, tspan, p1;
#     bcresid_prototype = (zeros(2), zeros(4))
# )
# sol_bvp = solve(prob, MIRK4(); dt = 0.05)
