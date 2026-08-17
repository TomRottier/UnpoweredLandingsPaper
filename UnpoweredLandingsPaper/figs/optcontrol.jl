export optcontrol

function optcontrol()
    f = Figure(size = (full_page_fig_width, 0.75full_page_fig_height), figure_padding = 5)
    cs = [:black; cs1[2:4]; cs1[1]]

    νtds = [1.0e-5, 0.1, 0.3, 0.5]
    extend = sqrt(1 / 0.12)
    xticks = ([0, extend], ["0", L"0.5 \nu_T^*"])

    ##### drag only
    pa = TrigAerodynamics(0.0, 0.03, 2.0)
    pc = ConstantAlpha(π / 2)
    p = Parameters(pa, pc)
    cb = ContinuousCallback((u, t, int) -> u[1] - extend, int -> terminate!(int))
    sols = map(enumerate(νtds)) do (i, νtd)
        γtd = asin(-max_CD(p) * νtd^2)
        return simulate(full_dynamics!, [νtd, γtd, 0.0, (i - 1)], p, (0.0, -1.0); callback = cb)
    end
    ax1a = phasespace_axis(
        [
            phasespace_regions(p; extend);
            mapreduce(vcat, enumerate(sols)) do (i, sol)
                trajectory_phasespace(sol; color = cs[i])
            end
        ]; extend, xticks, xlabel = "", xticklabelsvisible = false
    )

    ax1b = trajectory_axis(
        mapreduce(vcat, enumerate(sols)) do (i, sol)
            trajectory(sol, 1; τs = [0.0], color = cs[i])
        end;
        xlabel = "", xticklabelsvisible = false
    )
    gl1 = S.GridLayout([(1, 1, TopLeft()) => S.Label(text = "A", font = :bold, fontsize = 14pt, halign = 1.5, height = 12pt), (1, 1) => S.GridLayout([ax1a ax1b])])


    ###### lift only
    pa = TrigAerodynamics(1.5, 0.0, 0.0)
    pc = ConstantAlpha(π / 4)
    p = Parameters(pa, pc)
    cb = ContinuousCallback((u, t, int) -> u[2], nothing, int -> terminate!(int))
    sols = map(enumerate(νtds)) do (i, νtd)
        γtd = 0.0
        return simulate(full_dynamics!, [νtd, γtd, 0.0, (i - 1)], p, (0.0, -2.0); callback = cb)
    end
    # add trajectory for > 0 vertical climb speeds
    push!(sols, simulate(full_dynamics!, [0.7, π / 2, 0.0, -1.0], p, (0.0, -3.0); callback = cb))
    sol2 = simulate(full_dynamics!, [0.7, π / 2, 0.0, -1.0], Parameters(pa, ConstantAlpha(π / 2)), (0.0, 2.0); callback = ContinuousCallback((u, t, int) -> u[1] - 0.01, int -> terminate!(int)))

    ax2a = phasespace_axis(
        [
            phasespace_regions(p; extend);
            S.Lines(Point2f[sols[end].u[end][1:2], (0.0, π / 2)]; color = cs[end], linestyle = :dash);
            mapreduce(vcat, enumerate(sols)) do (i, sol)
                trajectory_phasespace(sol; color = cs[i])
            end
        ]; extend, xticks, xlabel = "", xticklabelsvisible = false
    )

    ax2b = trajectory_axis(
        [
            S.Lines(Point2f[sols[end][3:4, end], sol2[3:4, end]]; color = cs[end], linestyle = :dash);
            mapreduce(vcat, enumerate(sols)) do (i, sol)
                trajectory(sol, 1; τs = [0.0], color = cs[i])
            end
        ]; xlabel = "", xticklabelsvisible = false
    )
    gl2 = S.GridLayout([(1, 1, TopLeft()) => S.Label(text = "B", font = :bold, fontsize = 14pt, halign = 1.5, height = 12pt), (1, 1) => S.GridLayout([ax2a ax2b])])


    ###### lift and drag
    pa = TrigAerodynamics(1.5, 0.03, 2.0)
    pc = OptControl()
    p = Parameters(pa, pc)
    cb = ContinuousCallback((u, t, int) -> ν̇(u, int.p, t) + 0.01, int -> terminate!(int), nothing)
    sols = map(enumerate(νtds)) do (i, νtd)
        γtd = asin(-max_CD(p) * νtd^2) + 0.01
        return simulate(full_dynamics!, [νtd, γtd, 0.0, (i - 1)], p, (0.0, -3.0); callback = cb)
    end
    # add trajectory for > 0 vertical climb speeds
    push!(sols, simulate(full_dynamics!, [0.7, π / 2, 0.0, -1.0], p, (0.0, -3.0); callback = cb))
    sol2 = simulate(full_dynamics!, [0.7, π / 2, 0.0, -1.0], Parameters(pa, ConstantAlpha(π / 2)), (sols[end].t[end], 5.0); callback = ContinuousCallback((u, t, int) -> u[1] - 0.01, int -> terminate!(int)))

    # hypothetical dive trajectory - hermite spline for physical space traj, quadratic for phase space traj
    ν₁, γ₁, χ₁, ζ₁ = sols[1](0.0)
    χ̇₁, ζ̇₁ = sols[1](0.0, Val{1}; idxs = 3:4)
    χ₀, ζ₀, χ̇₀, ζ̇₀ = -5.5, sols[4][4, 1], 2.0, -3.0
    sol_dive_traj = map(0:0.01:1) do t
        x = (2t^3 - 3t^2 + 1) * χ₀ + (t^3 - 2t^2 + t) * χ̇₀ + (-2t^3 + 3t^2) * χ₁ + (t^3 - t^2) * χ̇₁
        z = (2t^3 - 3t^2 + 1) * ζ₀ + (t^3 - 2t^2 + t) * ζ̇₀ + (-2t^3 + 3t^2) * ζ₁ + (t^3 - t^2) * ζ̇₁
        Point2(x, z)
    end
    ν₀, γ₀ = sols[end - 1](0.0; idxs = 1:2)
    ν₂ = 0.5(ν₁ + ν₀)
    γ₂ = -1.0
    a, b, c = [ν₀^2 ν₀ 1; ν₁^2 ν₁ 1; ν₂^2 ν₂ 1] \ [γ₀, γ₁, γ₂]

    # colour white region
    undef_ν = range(sqrt(1 / max_CD(p)), sqrt(1 / min_CD(p)), 200)
    undef_γ = [acceleration_nullcline(ν, p) for ν in undef_ν]

    ax3a = phasespace_axis(
        [
            phasespace_regions(p; extend);
            S.Band(undef_ν, -π / 2, undef_γ, color = :white);
            S.Lines(undef_ν, undef_γ, color = cs3[2]);
            S.Lines(Point2f[sols[end].u[end][1:2], (0.0, π / 2)]; color = cs[end], linestyle = :dash);

            S.Lines(range(ν₀, ν₁, 100), ν -> evalpoly(ν, (c, b, a)); color = cs[end - 1], linestyle = (:dot, :dense));
            mapreduce(vcat, enumerate(sols)) do (i, sol)
                trajectory_phasespace(sol; color = cs[i])
            end
        ]; extend, xticks
    )

    ax3b = trajectory_axis(
        [
            S.Lines(sol_dive_traj; color = cs[end - 1], linestyle = (:dot, :dense));
            S.Scatter(χ₀, ζ₀; marker = bird, markersize = 30, color = cs[end - 1], rotation = -0.6);
            S.Lines(Point2f[sols[end][3:4, end], sol2[3:4, end]]; color = cs[end], linestyle = :dash);
            mapreduce(vcat, enumerate(sols)) do (i, sol)
                trajectory(sol, 1; τs = [0.0], color = cs[i])
            end
        ]; limits = (-7.0, 4.4, -2.4, 3.3), xticklabelspace = 12.0
    )

    gl3 = S.GridLayout([(1, 1, TopLeft()) => S.Label(text = "C", font = :bold, fontsize = 14pt, halign = 1.5, height = 12pt), (1, 1) => S.GridLayout([ax3a ax3b])])

    gl4 = S.GridLayout(
        [
            (1, 1, TopLeft()) => S.Label(text = "D", font = :bold, fontsize = 14pt, halign = 1.5, height = 12pt, tellheight = false, valign = 0.0),
            (1, 1) => timeseries_axis(
                [
                    mapreduce(vcat, enumerate(sols)) do (i, sol)
                        variable_timeseries(sol, 1; color = cs[i])
                    end
                    variable_timeseries(sol2, 1; color = cs[end], linestyle = :dash, decorate = false)
                ],
                L"speed $v / \hat{v}$"; yautolimitmargin = (0.0, 0.05), yticks = [0.0, 1.0, 2.0], limits = (nothing, nothing, 0, nothing)
            ),
            (1, 2) => timeseries_axis(
                [
                    mapreduce(vcat, enumerate(sols)) do (i, sol)
                        variable_timeseries(sol, 2; color = cs[i])
                    end
                    variable_timeseries(sol2, 2; color = cs[end], linestyle = :dash, decorate = false)
                ],
                L"flight path angle $\gamma$", yticks = ([0, π / 4, π / 2], ["0", L"\frac{π}{4}", L"\frac{π}{2}"])
            ),
            (1, 3) => timeseries_axis(
                [
                    S.HLines([π / 4, π / 2]; color = :gray, linestyle = :dash);
                    S.Text(
                        Point2f[(0.4, 0.9), (0.85, 1.4)];
                        text = [L"C_L = C_{L\text{max}}", L"C_D = C_{D\text{max}}"], fontsize = 9pt, align = (0.5, 0.5)
                    );
                    mapreduce(vcat, enumerate(sols)) do (i, sol)
                        alpha_timeseries(sol; color = cs[i])
                    end
                    alpha_timeseries(sol2; color = cs[end], linestyle = :dash, decorate = false)
                ], L"angle of attack $\alpha$", yticks = ([0, π / 4, π / 2], ["0", L"\frac{π}{4}", L"\frac{π}{2}"])
            ),
        ]
    )

    gl = S.GridLayout(
        [gl1, gl2, gl3, gl4];
        xaxislinks = [[ax1a, ax2a, ax3a], [ax1b, ax2b, ax3b]], yaxislinks = [[ax1a, ax2a, ax3a], [ax1b, ax2b, ax3b]]
    )

    plot(f[1, 1], gl)

    # align tick labels
    yspace = maximum(tight_yticklabel_spacing!, f.content[[8, 11]])
    [f.content[i].yticklabelspace = yspace for i in [8, 11]]
    f.content[8].xticklabelspace[] = f.content[9].xticklabelspace[]
    [f.content[i].yticklabelspace = 8.0 for i in [3, 6, 9]]

    map(i -> colgap!(f.layout.content[1].content.content[i].content.content[2].content, 1, 8), 1:3)

    return f

end
