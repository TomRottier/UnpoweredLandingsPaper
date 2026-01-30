using GlideLandingPaper, PointMass, OrdinaryDiffEqTsit5, CairoMakie
using ForwardDiff: gradient
using GlideLandingPaper: full_page_fig_height, full_page_fig_width, pt

∇α(u, p, t) = gradient(u -> PointMass.alpha(u, p, t), u)
∇α̇(u, p, t) = gradient(u -> ∇α(u, p, t)' * dynamics(u, p, t), u)
alpha(u, p, t) = u[5] - u[2]
∇γ̇(u, p, t) = gradient(u -> dynamics(u, p, t)[2], u)

function calculate_moment_arm(u, p, t)
    ν, γ, _, _, θ = u

    α = θ - γ
    cl = CL(α, p)
    cd = CD(α, p)
    ν̇ = -cd * ν^2 - sin(γ)
    γ̇ = cl * ν - cos(γ) / ν
    u̇ = [ν̇, γ̇]

    α̈ = ∇α̇([ν, γ], p, t)' * u̇
    γ̈ = ∇γ̇([ν, γ], p, t)' * u̇
    λ = (α̈ + γ̈) / (ν^2 * (cl * cos(α) + cd * sin(α)))

    return sign(λ) * min(10, abs(λ))
end

function full_dynamics(u, p, t)
    ν, γ, _, _, θ, θ̇ = u

    α = θ - γ
    cl = CL(α, p)
    cd = CD(α, p)
    λ = calculate_moment_arm(u, p, t)

    ν̇ = -cd * ν^2 - sin(γ)
    γ̇ = cl * ν - cos(γ) / ν
    χ̇ = ν * cos(γ)
    ζ̇ = ν * sin(γ)

    θ̈ = λ * ν^2 * (cl * cos(α) + cd * sin(α))

    return [ν̇, γ̇, χ̇, ζ̇, θ̇, θ̈]
end

function rb_comparison()
    f = Figure(size = (full_page_fig_width, 0.55full_page_fig_height))
    ax1 = Axis(f[1, 1][1, 1]; xlabel = L"range $x / \hat{l}$", ylabel = L"altitude $z / \hat{l}$", autolimitaspect = 1)
    ax2a = Axis(f[1, 2][1, 1]; xlabel = L"time $t / \hat{t}$", ylabel = L"angle of attack $\alpha$")
    ax2b = Axis(f[1, 2][2, 1]; xlabel = L"time $t / \hat{t}$", ylabel = L"moment arm $r / \hat{r}$")
    ax2c = Axis(f[1, 2][3, 1]; xlabel = L"time $t / \hat{t}$", ylabel = L"speed $v / \hat{v}$")

    ylims!(ax2a, 0, π / 2)
    ylims!(ax2b, -10, 10)
    linestyle = (:dash, :dense)

    νtds = [0.1, 0.3, 0.5]
    pa = TrigAerodynamics(1.5, 0.03, 2.0)
    pc = OptControl()
    p = Parameters(pa, pc)
    cb1 = ContinuousCallback((u, t, int) -> ν̇(u, int.p, t) + 0.02, int -> terminate!(int), nothing)

    for νtd in νtds
        γtd = -asin(max_CD(pa) * νtd^2)
        u₀_pm = [νtd, γtd, 0.0, 0.0]
        sol_pm = simulate(PointMass.full_dynamics!, u₀_pm, p, (0.0, -2.5); callback = cb1)

        θ₀ = sol_pm(0.0; idxs = 2) + PointMass.alpha(sol_pm(0.0), p, 0.0)
        θ̇₀ = ∇α(sol_pm(0.0; idxs = 1:2), p, 0.0)' * dynamics(sol_pm(0.0), p, 0.0) + dynamics(sol_pm(0.0), p, 0.0)[2]
        u₀_rb = [sol_pm(0.0); θ₀; θ̇₀]
        sol_rb = simulate(full_dynamics, u₀_rb, p, (0.0, 2.5); callback = cb1)

        lines!(ax1, sol_pm[3, :], sol_pm[4, :], label = L"\nu_{TD} = %$νtd")
        lines!(ax1, sol_rb[3, :], sol_rb[4, :]; linestyle)
        lines!(ax2a, sol_pm.t, t -> PointMass.alpha(sol_pm(t), p, t))
        l2 = lines!(ax2a, sol_rb.t, t -> alpha(sol_rb(t), p, t); linestyle)
        lines!(ax2b, sol_rb.t, t -> calculate_moment_arm(sol_rb(t), p, t); linestyle, color = l2.color[])
        lines!(ax2c, sol_pm.t, sol_pm[1, :])
        lines!(ax2c, sol_rb.t, sol_rb[1, :]; linestyle)

    end

    # labels
    Label(f[1, 1][1, 1, TopLeft()], "A"; font = :bold, fontsize = 14pt)
    Label(f[1, 2][1, 1, TopLeft()], "B"; font = :bold, fontsize = 14pt)
    Label(f[1, 2][2, 1, TopLeft()], "C"; font = :bold, fontsize = 14pt)
    Label(f[1, 2][3, 1, TopLeft()], "D"; font = :bold, fontsize = 14pt)

    rowgap!(f.layout.content[2].content, 0)
    yspace = maximum(tight_yticklabel_spacing!, [ax2a, ax2b, ax2c])
    [ax.yticklabelspace = yspace for ax in [ax2a, ax2b, ax2c]]

    return f
end
