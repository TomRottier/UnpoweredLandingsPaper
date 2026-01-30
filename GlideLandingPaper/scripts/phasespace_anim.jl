# animate how the phase space changes for a non-constant control strategy
using GlideLandingPaper, GLMakie, PointMass, OrdinaryDiffEqTsit5
using GlideLandingPaper: cs3, bird

pa = TrigAerodynamics(1.5, 0.03, 2.0)
pc = OptControl()
p = Parameters(pa, pc)
ν₁ = 0.2
γ₁ = -asin(max_CD(pa) * ν₁^2)
u₀ = [ν₁, γ₁, 0.0, 0.0]

function trajectory_gui(u₀, p; arrows = true)
    cb = ContinuousCallback((u, t, int) -> ν̇(u, int.p, t) + 0.1, int -> terminate!(int), nothing)
    sol = simulate(full_dynamics!, u₀, p, (0.0, -5.0); callback = cb)

    ### gui
    yulim = π / 2
    yllim = -π / 2
    extend = sqrt(1 / 0.03)
    n = 100
    f = Figure()
    ax1 = Axis(f[1:2, 1]; xlabel = L"speed $v / \hat{v}$", ylabel = L"flight path angle $\gamma$", limits = (0.0, extend, yllim, yulim))
    ax2 = Axis(
        f[1, 2]; xlabel = L"range $x / \hat{d}$", ylabel = L"altitude $z / \hat{d}$",
        autolimitaspect = 1, xautolimitmargin = (0.2, 0.2)
    )
    ax3 = Axis(f[2, 2], xlabel = L"time $t / \hat{t}$", ylabel = L"\alpha")

    sg = SliderGrid(
        f[3, :],
        (label = L"time $t / \hat{t}$", range = sol.t[begin]:0.001:(sol.t[end]), startvalue = sol.t[begin])
    )
    τ = lift(sg.sliders[1].value) do t
        t
    end

    ν = @lift sol($τ)[1]
    γ = @lift sol($τ)[2]
    α = @lift alpha([$ν, $γ], p, $τ)
    θ = @lift $α + $γ
    pctemp = @lift ConstantAlpha($α)
    ptemp = @lift Parameters(pa, $pctemp)

    if arrows
        us = Point2f[(ν, γ) for ν in range(0.0, extend; length = 30), γ in range(-π / 2, π / 2; length = 30)]
        ups = @lift Point2f[dynamics(u, $ptemp, nothing) for u in us]
        arrows2d!(ax1, us, ups; lengthscale = 0.1, color = :black, strokemask = 0.0)
    end

    # acceleration nullcline
    νs = @lift acceleration_nullcline($ptemp; n)[1]
    γs = @lift acceleration_nullcline($ptemp; n)[2]
    lower = @lift Point2.($νs, $γs)
    upper = @lift Point2.($νs, yulim)
    band!(ax1, lower, upper, color = cs3[1], alpha = 0.5)

    lower = @lift Point2f[($νs[end], yllim), (extend, yllim)]
    upper = @lift Point2f[($νs[end], yulim), (extend, yulim)]
    band!(ax1, lower, upper; color = cs3[1], alpha = 0.5)
    lines!(ax1, νs, γs, color = cs3[2])

    lower = @lift Point2.($νs, -π / 2)
    upper = @lift Point2.($νs, $γs)
    band!(ax1, lower, upper; color = cs3[2], alpha = 0.5)

    # flightpath nullcline
    νs = @lift flightpath_nullcline($ptemp; n)[1]
    γs = @lift flightpath_nullcline($ptemp; n)[2]
    xs = @lift [$νs; reverse($νs)]
    ys = @lift [$γs; -reverse($γs)]
    lines!(ax1, xs, ys, color = cs3[3])

    lower = @lift Point2.($νs[1:(end ÷ 2)], -1 * $γs[1:(end ÷ 2)])
    upper = @lift Point2.($νs[1:(end ÷ 2)], $γs[1:(end ÷ 2)])
    band!(ax1, lower, upper, color = cs3[3], alpha = 0.5)


    # trajectory
    pts = @lift Point2(sol($τ; idxs = 3:4))
    lines!(ax2, [Point2(χ, ζ) for (_, _, χ, ζ) in sol.u]; color = :gray)
    scatter!(ax2, pts; marker = bird, markersize = 50, rotation = θ)

    pts = @lift Point2($ν, $γ)
    lines!(ax1, [Point2(ν, γ) for (ν, γ) in sol.u]; color = :gray)
    scatter!(ax1, pts)

    lines!(ax3, sol.t, t -> alpha(sol(t), p, t))
    scatter!(ax3, τ, α)

    return f
end


### makie video file
using CairoMakie
CairoMakie.activate!()

function create_video(; ft = ".mp4")
    pa = TrigAerodynamics(1.5, 0.03, 2.0)
    pc = OptControl()
    p = Parameters(pa, pc)
    ν₁ = 0.2
    γ₁ = -asin(max_CD(pa) * ν₁^2)
    u₀ = [ν₁, γ₁, 0.0, 0.0]
    cb = ContinuousCallback((u, t, int) -> ν̇(u, int.p, t) + 0.1, int -> terminate!(int), nothing)
    sol = simulate(full_dynamics!, u₀, p, (0.0, -5.0); callback = cb)

    yulim = π / 2
    yllim = -π / 2
    extend = sqrt(1 / 0.03)

    f = Figure()
    ax1 = Axis(f[1:2, 1]; xlabel = L"speed $v / \hat{v}$", ylabel = L"flight path angle $\gamma$", limits = (0.0, extend, yllim, yulim))
    ax2 = Axis(
        f[1, 2]; xlabel = L"range $x / \hat{d}$", ylabel = L"altitude $z / \hat{d}$",
        autolimitaspect = 1, xautolimitmargin = (0.2, 0.2)
    )
    ax3 = Axis(f[2, 2], xlabel = L"time $t / \hat{t}$", ylabel = L"\alpha", yticks = ([0, π / 4, π / 2], ["0", "π/4", "π/2"]))

    lines!(ax1, [Point2(ν, γ) for (ν, γ) in sol.u]; color = :gray)
    lines!(ax2, [Point2(χ, ζ) for (_, _, χ, ζ) in sol.u]; color = :gray)
    lines!(ax3, sol.t, t -> alpha(sol(t), p, t); color = :gray)

    N = 200
    τs = range(sol.t[begin], sol.t[end]; length = N)

    return record(f, "figs//phasespace" * ft, τs, framerate = 24) do τ
        map(ax -> [delete!(ax, plt) for plt in ax.scene.plots[2:end]], [ax1, ax2, ax3])
        ν = sol(τ)[1]
        γ = sol(τ)[2]
        α = alpha([ν, γ], p, τ)
        θ = α + γ
        pctemp = ConstantAlpha(α)
        ptemp = Parameters(pa, pctemp)

        # acceleration nullcline
        νs = acceleration_nullcline(ptemp; n = 100)[1]
        γs = acceleration_nullcline(ptemp; n = 100)[2]
        lower = Point2.(νs, γs)
        upper = Point2.(νs, yulim)
        band!(ax1, lower, upper, color = cs3[1], alpha = 0.5)

        lower = Point2f[(νs[end], yllim), (extend, yllim)]
        upper = Point2f[(νs[end], yulim), (extend, yulim)]
        band!(ax1, lower, upper; color = cs3[1], alpha = 0.5)
        lines!(ax1, νs, γs, color = cs3[2])

        lower = Point2.(νs, -π / 2)
        upper = Point2.(νs, γs)
        band!(ax1, lower, upper; color = cs3[2], alpha = 0.5)

        # flightpath nullcline
        νs = flightpath_nullcline(ptemp; n = 100)[1]
        γs = flightpath_nullcline(ptemp; n = 100)[2]
        xs = [νs; reverse(νs)]
        ys = [γs; -reverse(γs)]
        lines!(ax1, xs, ys, color = cs3[3])

        lower = Point2.(νs[1:(end ÷ 2)], -1 * γs[1:(end ÷ 2)])
        upper = Point2.(νs[1:(end ÷ 2)], γs[1:(end ÷ 2)])
        band!(ax1, lower, upper, color = cs3[3], alpha = 0.5)

        # trajectory
        pts_p = Point2(sol(τ; idxs = 3:4))
        pts_v = Point2(ν, γ)
        arrows2d!(ax2, pts_p, pts_v; lengthscale = 0.25)
        scatter!(ax1, pts_v)
        scatter!(ax2, pts_p; marker = bird, markersize = 50, rotation = θ)
        scatter!(ax3, τ, α)
    end
end
