export freebodydiagram

function freebodydiagram()
    f = Figure(size=(full_page_fig_width, 0.75full_page_fig_height), figure_padding=(5, 8, 5, 5))

    α = π / 4
    pa = TrigAerodynamics(1.5, 0.03, 2.0)
    p = Parameters(pa, ConstantAlpha(α))
    ν = 1.2
    γ = 0.4
    cl = CL(α, p)
    cd = CD(α, p)
    sf = 0.27

    cb1 = ContinuousCallback((u, t, int) -> u[2] + 0.4, int -> terminate!(int))
    cb2 = ContinuousCallback((u, t, int) -> ν̇(u, int.p, t), int -> terminate!(int), nothing)
    sol1 = simulate(full_dynamics!, [ν, γ, 0.0, 0.0], p, (0.0, -5.0); callback=cb1)
    sol2 = simulate(full_dynamics!, sol1.u[1], p, (0.0, 6.0); callback=cb2)
    sol3 = simulate(full_dynamics!, sol2.u[end], p, (0.0, 6.0))

    τfbd = sol1.t[end]
    τstall = solve(NonlinearProblem((t, _) -> γ̇(sol2(t), p, t), 0.5)).u
    τterm = sol2.t[end]
    νstall, γstall = sol2(τstall; idxs=1:2)
    νxstall, νzstall = νstall * cos(γstall), νstall * sin(γstall)
    Lstall = cl * νstall^2
    Lx, Lz = -Lstall * sin(γstall), Lstall * cos(γstall)
    νterm, γterm = sol2(τterm; idxs=1:2)
    νxterm, νzterm = νterm * cos(γterm), νterm * sin(γterm)
    Dterm = cd * νterm^2
    Dx, Dz = -Dterm * cos(γterm), -Dterm * sin(γterm)

    d = 0.25 # axis limits around pt for fbd

    idx1 = findfirst(≥(-0.5d), sol2[3, :])
    idx2 = findfirst(≥(0.5d), sol2[3, :])
    traj = [sol2(t; idxs=3:4) for t in sol2.t[idx1:idx2]]
    ax1 = diagram_axis(free_body_diagram(ν, γ, p, traj; sf, pt=[0.57, 0.45], gpt=[0.1, 0.8], upt=[0.95, 0.9], birdsize=60))


    ax2 = trajectory_axis(
        [
            trajectory(
                sol2, 4;
                τs=[0.05, abs(τfbd)], forces=false, weight=false, vel=true, sf=0.17, decorate_end=false, markersize=50
            )...,
            trajectory(
                sol3, 2;
                τs=[0.5, 1.3], forces=false, weight=false, vel=true, sf=0.17, decorate_start=false, markersize=50, linestyle=(:dot, :dense)
            )...,
            S.Poly(Rect2f((-d / 2, -d / 2), (d, d)), color=:transparent, strokewidth=1),
            S.Scatter(Point2(sol2(τstall; idxs=3:4)); marker=:rect, color=cs1[3]),
            S.Scatter(Point2(sol2(τterm; idxs=3:4)); marker=:rect, color=cs1[2]),
            S.Text(Point2(sol2(τstall; idxs=3:4)); text="D", align=(:right, :bottom), font=:bold),
            S.Text(Point2(sol2(τterm; idxs=3:4)); text="E", align=(:right, :bottom), font=:bold),
            S.Text(Point2(0, 0); text="A", align=(:right, :bottom), font=:bold, offset=(0, 20)),
        ];
        xticksvisible=false, xticklabelsvisible=false,
        yticksvisible=false, yticklabelsvisible=false, yticklabelspace=15.0,
        alignmode=Inside()
    )

    ax3 = phasespace_axis(
        [
            phasespace_regions(p; extend=sqrt(1 / 0.15))...,
            trajectory_phasespace(sol2, decorate_end=false)...,
            trajectory_phasespace(sol3; decorate_start=false, linestyle=(:dot, :dense))...,
            S.Scatter(Point2(sol2(τstall; idxs=1:2)); marker=:rect, color=cs1[3]),
            S.Scatter(Point2(sol2(τterm; idxs=1:2)); marker=:rect, color=cs1[2]),
            S.Scatter(Point2(ν, γ); marker=:rect, color=:grey),
            S.Text(
                Point2f[(0.65, 0.75), (0.1, 0.7), (0.03, 0.25), (0.25, 0.05)];
                text=[L"\dot{\nu} < 0, \, \dot{\gamma} > 0", L"\dot{\nu} < 0, \, \dot{\gamma} < 0", L"\dot{\nu} > 0, \, \dot{\gamma} < 0", L"\dot{\nu} > 0, \, \dot{\gamma} > 0"], space=:relative
            ),
            S.Text(Point2(sol2(τstall; idxs=1:2)); text="D", align=(:left, :bottom), font=:bold),
            S.Text(Point2(sol2(τterm; idxs=1:2)); text="E", align=(:left, :bottom), font=:bold),
            S.Text(Point2(ν, γ); text="A", align=(:left, :bottom), font=:bold),
            S.Text(Point2(νstall, γstall); text=L"\nu_{\text{stall}}", offset=(-10, 23), align=(:center, :baseline), color=:black, fontsize=14pt),
            S.Text(Point2(sqrt(1 / cd), -π / 2); text=L"\nu_{\text{term}}", offset=(15, 20), align=(:center, :baseline), color=:black, fontsize=14pt), # to_color((1.0, 192 / 255, 114 / 255))
        ];
        extend=sqrt(1 / 0.15), xticklabelspace=5.0, yticklabelspace=15.0, alignmode=Inside()
    )

    sf = 0.54
    ax4a = diagram_axis(
        [
            stall_speed_diagram(νstall, γstall, p; sf, x=0.3, y=0.6, arrow_alpha=0.3, textx=0.7, texty=0.15, forces=false, narrows=15, markersize=60);
        ]
    )

    ax4b = diagram_axis(
        [
            term_speed_diagram(νterm, γterm, p; sf, x=0.3, y=0.6, arrow_alpha=0.3, textx=0.7, texty=0.85, forces=false, narrows=15, markersize=60);
        ]
    )

    ν, γ = 1.0, 0.0
    kwargs = (x=0.25, y=0.65, sf=0.4, forces=false, fontsize=0, arrow_alpha=0.1, markersize=60, fbd_arrows=false, narrows=15)
    ax4c = S.GridLayout(
        [
        (1, 1) => diagram_axis(
            [
                stall_speed_diagram(ν, γ, Parameters(pa, ConstantAlpha(0.4)); kwargs...);
                term_speed_diagram(ν, γ, Parameters(pa, ConstantAlpha(0.4)); kwargs...)
            ]
        ),
        (1, 1) => diagram_axis(
            [
                stall_speed_diagram(ν, γ, Parameters(pa, ConstantAlpha(0.5)); kwargs...);
                term_speed_diagram(ν, γ, Parameters(pa, ConstantAlpha(0.5)); kwargs...)
            ]
        ),
        (1, 1) => diagram_axis(
            [
                stall_speed_diagram(ν, γ, Parameters(pa, ConstantAlpha(π / 4)); kwargs...);
                term_speed_diagram(ν, γ, Parameters(pa, ConstantAlpha(π / 4)); kwargs...);
                S.Arrows2D(Point2(0.45, 0.75), Point2(0.35, 0.0); color=cs3[3], space=:relative);
                S.Text(Point2(0.8, 0.75); text=L"\downarrow C_L", color=cs3[3], align=(:left, :center), space=:relative);
                S.Arrows2D(Point2(0.35, 0.35), Point2(0.35, -0.2); color=cs3[2], space=:relative);
                S.Text(Point2(0.7, 0.15); text=L"\downarrow C_D", color=cs1[2], align=(:left, :center), space=:relative)
            ]
        ),
    ]
    )

    CairoMakie.plot(
        f[1, 1],
        S.GridLayout(
            [
            (1, 1) => S.GridLayout(
                [
                (1, 1) => ax1,
                (1, 1, TopLeft()) => S.Label(; text="A", fontsize=14pt, font=:bold),
            ]
            ),
            (1, 2) => S.GridLayout(
                [
                (1, 1) => ax2,
                (1, 1, TopLeft()) => S.Label(; text="B", fontsize=14pt, font=:bold),
            ]
            ),
            (2, :) => S.GridLayout(
                [
                (1, 1) => ax3,
                (1, 1, TopLeft()) => S.Label(; text="C", fontsize=14pt, font=:bold, padding=(0, 0, 0, 0)),
            ]
            ),
            (3, :) => S.GridLayout(
                [
                (1, 1) => ax4a,
                (1, 1, TopLeft()) => S.Label(; text="D", fontsize=14pt, font=:bold),
                (1, 2) => ax4b,
                (1, 2, TopLeft()) => S.Label(; text="E", fontsize=14pt, font=:bold),
                (1, 3) => ax4c,
                (1, 3, TopLeft()) => S.Label(; text="F", fontsize=14pt, font=:bold),
            ]
            ),
        ]
        )
    )


    rowgap!(f.layout.content[1].content, 0)
    colsize!(f.layout.content[1].content, 1, Aspect(1, 1.0))
    rowsize!(f.layout.content[1].content.content[4].content, 1, Aspect(1, 1.0))
    limits!(f.content[3], nothing, 1.3, -0.2, nothing)
    linkaxes!(f.content[[7; 9]]...)
    linkaxes!(f.content[11:13]...)
    resize_to_layout!(f)


    return f
end


# # connect fbd to trajectory
# _ax1 = f.content[1] # the actual Axis object
# _ax2 = f.content[2]
# _ax3 = f.content[3]
# pt11 = viewport(_ax1.scene)[].origin + [viewport(_ax1.scene)[].widths[1], 0]
# pt21 = pt11 + [0, viewport(_ax1.scene)[].widths[1]]
# pt12 = Point2f(-0.5d, -0.5d)
# pt22 = Point2f(-0.5d, 0.5d)
# lines!(f.scene, Point2f[pt11, Makie.shift_project(_ax2.scene, pt12)];
#     color=:black, linestyle=(:dot, :dense))
# lines!(f.scene, Point2f[pt21, Makie.shift_project(_ax2.scene, pt22)];
#     color=:black, linestyle=(:dot, :dense))

# # connect fbd to phase space
# pt31 = viewport(_ax1.scene)[].origin
# pt3 = Point2f(ν, γ)
# l1 = lines!(f.scene, Point2f[pt31, Makie.shift_project(_ax3.scene, pt3)];
#     color=:black, linestyle=(:dot, :dense))
# l2 = lines!(f.scene, Point2f[pt11, Makie.shift_project(_ax3.scene, pt3)];
#     color=:black, linestyle=(:dot, :dense))
# translate!(l1, (0, 0, 10000))
# translate!(l2, (0, 0, 10000))
