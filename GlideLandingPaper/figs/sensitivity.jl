export sensitivity

"""
    get_initial_speed(νtd, p)

Simulates trajectory backwards from touchdown at speed `νtd` > 0 and `γtd` required to give ν̇ = 0 to a horizontal flight path for a given set of parameters `p`.
"""
function get_initial_speed(νtd, p)
    cbs = CallbackSet(
        ContinuousCallback((u, t, int) -> u[2], nothing, int -> terminate!(int)),
        ContinuousCallback((u, t, int) -> ν̇(u, int.p, t) + 0.05, int -> terminate!(int), nothing)
    )
    γtd = acceleration_nullcline(νtd, p)
    sol = simulate(dynamics!, [νtd, γtd], p, (0.0, -5.0); callback=cbs, save_everystep=false)
    return sol.u[1][1]
end

"""
    energy_height(ν, pa)

Calculate the initial height needed to reach an initial speed assuming diving vertically downwards at minimum CD.

Analytical solution to given differential equation `ν̇ = -CD₀ * ν² + 1` is given by:

    `ν = tanh(√CD₀ * τ) / √CD₀`
"""
function energy_height(ν, pa)
    CD₀ = min_CD(pa)
    ν ≥ 1 / √CD₀ && return 1000.0
    if CD₀ > 0
        τ = atanh(ν * √CD₀) / √CD₀ # time to reach ν
        Eloss = -(tanh(√CD₀ * τ)^2 - 2log(cosh(√CD₀ * τ))) / 2CD₀ # energy lost to drag
    else
        τ = ν # for small x tanh(x) ≈ x - x³/3 ∴ ν = τ
        Eloss = 0.0
    end

    return 0.5ν^2 + Eloss
end

function sensitivity()
    f = Figure(size=(full_page_fig_width, 0.75full_page_fig_height))

    νtds = [0.1, 0.4]
    νtermmin = sqrt(1 / 0.03)
    N = 51
    CLmaxs = range(0.1, 5.0, N)
    CDmaxs = range(0.1, 5.0, N)
    pc = OptControl()
    cbs = CallbackSet(
        ContinuousCallback((u, t, int) -> u[1] - νtermmin, int -> terminate!(int)),
        ContinuousCallback((u, t, int) -> ν̇(u, int.p, t) + 0.1, int -> terminate!(int), nothing)
    )

    pas = [
        TrigAerodynamics(1.5, 0.03, 2.0),
        TrigAerodynamics(3.0, 0.03, 2.0),
        TrigAerodynamics(1.5, 0.03, 4.0),
        TrigAerodynamics(3.0, 0.03, 4.0),
    ]
    cols = [to_color(:grey), cs1[1:3]...]
    levels = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, νtermmin]

    axs = map(νtds) do νtd
        sols = Matrix{ODESolution}(undef, N, N)
        Threads.@threads for (idx, (CLmax, CDmax)) in collect(enumerate(Iterators.product(CLmaxs, CDmaxs)))
            pa = TrigAerodynamics(CLmax, 0.03, CDmax)
            p = Parameters(pa, pc)
            γtd = asin(-max_CD(pa) * νtd^2) + 0.1
            u1 = [νtd, γtd, 0.0, 0.0]
            sols[idx] = simulate(full_dynamics!, u1, p, (0.0, -5.0); callback=cbs, save_on=false)
        end

        touchdown_speeds = map(sols) do sol
            sol.u[1][1]
        end

        selected_sols = map(pas) do pa
            p = Parameters(pa, OptControl())
            γtd = asin(-max_CD(pa) * νtd^2) + 0.1
            u1 = [νtd, γtd, 0.0, 0.0]
            simulate(full_dynamics!, u1, p, (0.0, -5.0); callback=cbs)
        end

        ax1 = S.Axis(
            plots=[
                S.Contourf(CLmaxs, CDmaxs, touchdown_speeds; colormap=Reverse(:RdBu), levels, extendhigh=to_colormap(:RdBu)[1]);
                S.Scatter([Point2f(p.CLmax, p.CDmax) for p in pas]; color=cols, marker=:rect, markersize=14)
            ]; xlabel=L"C_{L\text{max}}", ylabel=L"C_{D\text{max}}", autolimitaspect=1, title="touchdown speed $νtd", xautolimitmargin=(0.0, 0.0), yautolimitmargin=(0.0, 0.0), titlesize=12pt, titlefont=:regular
        )

        ax2 = trajectory_axis(
            reduce(vcat, trajectory(sol, 2; τs=[0.0], color=cols[i], markersize=45) for (i, sol) in enumerate(selected_sols))
        )

        return [ax1, ax2]
    end

    gl = S.GridLayout(vec([(j, i) => axs[i][j] for i in 1:2, j in 1:2]))

    plot(f[1, 1], gl)
    tightlimits!(f.content[1])
    tightlimits!(f.content[2])
    Colorbar(f[0, :], f.content[1].scene.plots[1]; ticks=([0, 1, 2, 3, 4, 5, 5.7], ["0", "1", "2", "3", "4", "5", L"\nu_{\text{max}}"]), vertical=false, tellheight=true, label=L"initial speed $v / \hat{v}$", labelsize=12pt)

    rowsize!(f.layout.content[1].content, 1, Aspect(1, 1.0))

    return f
end


#### effect of varying CD0
function cd0_sensitivity()
    cd0s = [0.03, 0.06, 0.015]
    νtd = 0.1
    νinits = map(cd0s) do cd0
        pa = TrigAerodynamics(1.5, cd0, 2.0)
        p = Parameters(pa, OptControl())
        return get_initial_speed(νtd, p)
    end

    # % change in νinit for each CD0
    return Dict(
        cd0 => νinit / νinits[1] for (cd0, νinit) in zip(cd0s, νinits)
    )
end
