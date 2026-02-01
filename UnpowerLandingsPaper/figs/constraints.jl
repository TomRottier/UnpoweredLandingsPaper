export constraints

function constraints()
    f = Figure(size = (full_page_fig_width, 0.67full_page_fig_height), figure_padding = 5)

    # generate solutions
    pa = TrigAerodynamics(1.5, 0.03, 2.0)
    ν₀ = 1.3
    γ₀ = -0.18
    u₀ = [ν₀, γ₀, 0.0, 0.0]
    νh = sqrt(1 / max_CL(pa))
    νv = sqrt(1 / max_CD(pa))

    init = nothing
    sol_min = solve_ocp(
        dc_min_speed_ground(u₀, pa; Tmin = 3.0, Tmax = 3.0); init, solver = :madnlp, display = false
    )
    νmin, γmin = sol_min.u[end][1:2]
    Tsol_min = sol_min.t[end]

    γs1 = range(0.0, γmin + 0.1, 6)
    sols1 = map(γs1) do γ
        solve_ocp(
            dc_min_speed_constrain_td_angle(u₀, pa, γ; Tmin = Tsol_min, Tmax = Tsol_min); init = sol_min.sol, solver = :madnlp, display = false
        )
    end

    γs2 = range(γmin - 0.1, -π / 2, 8)
    sols2 = map(γs2) do γ
        sol = solve_ocp(
            dc_min_speed_constrain_td_angle(u₀, pa, γ; Tmax = 3.0); init = sol_min.sol, solver = :madnlp, display = false
        )
        OptimalControl.successful(sol.sol) || @warn "unsuccessful at γ: $γ"
        sol
    end
    _sols3 = simulate(full_dynamics!, u₀, Parameters(pa, ConstantAlpha(π / 2)), (0.0, 3.0))
    sols3 = simulate(full_dynamics!, [u₀[1:3]..., -_sols3[4, end]], Parameters(pa, ConstantAlpha(π / 2)), (0.0, 3.0))

    γs = [γs1; γmin; γs2; -π / 2]
    sols = [sols1; sol_min; sols2; sols3]
    colorrange = (-π / 2, 0.0)
    νtds = map(sols) do sol
        T = sol.t[end]
        return sol(T)[1]
    end
    νmin, minidx = findmin(νtds)
    γmin = sols[minidx].u[end][2]
    special_sols = [1, minidx, length(sols)]


    ### trajectories
    ax1a = trajectory_axis(
        [
            S.HLines(0.0; color = :gray, linewidth = 1);
            mapreduce(vcat, enumerate(sols)) do (i, sol)
                pts = [Point2f(sol(t)[3:4]) for t in sol.t]
                return S.Lines(pts; alpha = 0.2, colorrange, color = γs[i])
            end
            mapreduce(vcat, special_sols) do i
                sol = sols[i]
                pts = [Point2f(sol(t)[3:4]) for t in sol.t]
                νtd, γtd, χtd, ζtd = sol.u[end]
                T = sol.t[end]
                l = S.Lines(pts; colorrange, color = γs[i])
                α = sol isa ODESolution ? alpha([νtd, γtd], sol.prob.p, T) : sol.cntrlf(T)
                s = S.Scatter(Point2f(χtd, ζtd); marker = bird, colorrange, color = γs[i], rotation = sol(T)[2] + α, markersize = 50)
                a = velocity_vector(Point2f(χtd, ζtd), νtd, γtd; colorrange, color = γs[i], sf = 1.0)
                return [l; a; s]
            end
        ]; xautolimitmargin = (0.0, 0.01), limits = (0.0, nothing, nothing, nothing), xticklabelspace = 10.0
    )

    ax2a = timeseries_axis(
        [
            S.HLines([sqrt(1 / max_CL(pa)), sqrt(1 / max_CD(pa))]; colorrange, color = γs[[1, end]], linestyle = :dash);
            S.Text(
                Point2f[(0.4, νh), (1.0, νv)];
                text = [L"\sqrt{1/C_{L\text{max}}}", L"\sqrt{1/C_{D\text{max}}}"], colorrange, color = γs[[1, end]], fontsize = 10pt, offset = (0, 3)
            );
            map(special_sols) do i
                sol = sols[i]
                pts = [Point2f(t, u[1]) for (u, t) in zip(sol.u, sol.t)]
                S.Lines(pts; colorrange, color = γs[i])
            end
        ],
        L"speed $v / \hat{v}$"; limits = (0.0, 3.0, nothing, nothing)
    )
    ax2b = timeseries_axis(
        [
            S.HLines([0.0, -π / 2]; colorrange, color = γs[[1, end]], linestyle = :dash);
            map(special_sols) do i
                sol = sols[i]
                pts = [Point2f(t, u[2]) for (u, t) in zip(sol.u, sol.t)]
                S.Lines(pts; colorrange, color = γs[i])
            end
        ],
        L"flight path angle $\gamma$"; limits = (0.0, 3.0, nothing, nothing), yticks = ([-π / 2, -π / 4, 0.0], [L"-$\frac{π}{2}$", L"-$\frac{π}{4}$", "0"])
    )
    ax2c = timeseries_axis(
        [
            S.HLines([π / 4, π / 2]; colorrange, color = γs[[1, end]], linestyle = :dash);
            map(special_sols) do i
                sol = sols[i]
                _f = sol isa ODESolution ? t -> alpha(sol(t), sol.prob.p, t) : sol.cntrlf
                pts = [Point2f(t, _f(t)) for (u, t) in zip(sol.u, sol.t)]
                S.Lines(pts; colorrange, color = γs[i])
            end
            S.Text([Point2(0.7, 0.65), Point2(1.25, 1.4)], text = [L"C_L = C_{L\text{max}}", L"C_D = C_{D\text{max}}"], fontsize = 11pt, align = (:center, :baseline), offset = (0, 0))
        ],
        L"angle of attack $\alpha$"; limits = (0, 3.0, 0, π / 2 + 0.05), yticks = ([0.0, π / 4, π / 2], ["0", L"\frac{π}{4}", L"\frac{π}{2}"])
    )

    gl = S.GridLayout(
        (0, 1:2) => S.Colorbar(; limits = (-π / 2, 0.0), colormap = cs2, vertical = false, label = "touchdown angle", ticks = ([-π / 2, -π / 4, 0.0], ["-π/2", "-π/4", "0"])),
        (1, 1:2) => S.GridLayout(
            [
                (1, 1) => ax1a,
                (1, 1, TopLeft()) => S.Label(; text = "A", fontsize = 14pt, font = :bold),
            ]
        ),
        (2, 1:2) => S.GridLayout(
            [
                (1, 1) => ax2a,
                (1, 1, TopLeft()) => S.Label(; text = "Ci", fontsize = 14pt, font = :bold),
                (1, 2) => ax2b,
                (1, 2, TopLeft()) => S.Label(; text = "Cii", fontsize = 14pt, font = :bold),
                (1, 3) => ax2c,
                (1, 3, TopLeft()) => S.Label(; text = "Ciii", fontsize = 14pt, font = :bold),
            ]
        )
    )

    plot(f[1, 1], gl)


    # plot old fashioned way as PolarAxis doesn't appear to work with SpecApi
    pax = PolarAxis(
        f[1, 1][1, 1:2];
        thetalimits = (-π / 2, 0.0), rlimits = (0, 1.0), rgridvisible = true, thetagridvisible = false, thetaticks = ([0, -π / 4, -π / 2, γmin], ["0", "-π/4", "-π/2", L"\gamma_{\text{min}}"]), rticks = [0.5, 1.0], spinevisible = false, rticklabelsvisible = true, rticklabelpad = 5, clip = false, height = Relative(0.9), width = Relative(0.9), valign = 0.3, halign = 3.5, title = L"touchdown speed $v / \hat{v}$", titlesize = 15, titlegap = 0, backgroundcolor = :white
    )
    lines!(pax, Point2f[(0.0, 0.0), (0.0, 1.0), (NaN, NaN), (0.0, 0.0), (-π / 4, 1.0), (NaN, NaN), (0.0, 0.0), (-π / 2, 1.0)]; color = :gray, linewidth = 1)
    lines!(pax, Point2f[(0.0, 0.0), (γmin, 1.0)]; color = :gray, linestyle = :dash, linewidth = 1)
    lines!(pax, -π / 2 .. 0.0, γ -> sqrt(-sin(γ) / max_CD(pa)); color = cs3[2], alpha = 0.6, label = L"\nu_{\text{term},\, C_{D\text{max}}}")
    lines!(pax, -π / 2 .. 0.0, γ -> sqrt(cos(γ) / max_CL(pa)); color = cs3[3], alpha = 0.6, label = L"\nu_{\text{stall},\, C_{L\text{max}}}")
    arrows2d!(pax, fill(Point2f(0.0, 0.0), length(special_sols)), [Point2f(γs[i], 1.07νtds[i]) for i in special_sols]; colorrange, color = γs[special_sols], tip = arrowhead_vel, shaftwidth = 3, tipwidth = 10, tiplength = 10)
    scatter!(pax, γs, νtds; colorrange, color = γs) # linecap = :square)
    θs = map(sols[special_sols]) do sol
        _f = sol isa ODESolution ? t -> alpha(sol(t), sol.prob.p, t) : sol.cntrlf
        return sol.u[end][2] + _f(sol.t[end])
    end
    scatter!(pax, fill(Point2(0, 0), length(special_sols)); marker = bird, markersize = 65, colorrange, color = γs[special_sols], rotation = θs, alpha = 1.0)
    text!(
        pax, Point2f[(-0.3, 1.15sqrt(cos(-0.3) / max_CL(pa))), (-1.1, 1.25sqrt(-sin(-1.1) / max_CD(pa)))];
        text = [L"\nu_{\text{stall}}", L"\nu_{\text{term}}"], align = (:center, :baseline), color = [cs3[3], to_color((238 / 255, 174 / 255, 24 / 255))], fontsize = 12pt
    )
    Label(f[1, 1][1, 1:2]; text = "B", fontsize = 14pt, font = :bold, valign = 0.95, halign = 0.6)


    yspace = maximum(tight_yticklabel_spacing!, f.content[[2, 4]])
    [f.content[i].yticklabelspace = yspace for i in [2, 4]]
    rowsize!(f.layout.content[1].content, 1, Relative(0.5))
    rowgap!(f.layout.content[1].content, 5)

    return f
end
