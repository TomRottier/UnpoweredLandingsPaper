export dynamics_phasespace

function dynamics_phasespace()
    pas = [TrigAerodynamics(1.5, 0.0, 0.0), TrigAerodynamics(1.5, 0.03, 2.0)]
    pcs = [ConstantAlpha(0.0), ConstantAlpha(0.121), ConstantAlpha(π / 4), ConstantAlpha(π / 2)]
    u₀s = [
        [(0.15, 0.0), (0.9, -0.4), (1.3, 0.5)],
        [(0.15, 0.0), (0.9, -0.4), (1.3, 0.5)],
        [(0.15, 0.0), (0.9, -0.4), (1.3, 0.5)],
        [(0.15, 0.0), (0.9, -0.4), (1.3, 0.5)],
    ]
    axlabels1 = ["A", "B", "C", "D"]
    axlabels2 = ["i", "ii"]
    gls = Pair{Tuple{Int, Int}, Makie.GridLayoutSpec}[]

    f = Figure(size = (full_page_fig_width, full_page_fig_height - 2inch), figure_padding = (5, 5, 5, 16pt))
    for (i, pc) in enumerate(pcs)
        for (j, pa) in enumerate(pas)
            p = Parameters(pa, pc)
            α = alpha(nothing, p, nothing)
            cl, cd = CL(α, p), CD(α, p)
            νstar = 1 / fourthroot(cl^2 + cd^2)
            γstar = atan(-cd / cl)
            extend = sqrt(1 / 0.03)
            xlabel = i == length(pcs) ? L"speed $v / \hat{v}$" : ""
            ylabel = j == 1 ? L"flight path angle $\gamma$" : ""
            xticks = i == length(pcs) ? ([0.0, extend], ["0", L"\nu_T^*"]) : ([0, extend], ["", ""])
            yticks = j == 1 ? ([-π / 2, 0.0, π / 2], [L"-$\frac{π}{2}$", "0", L"\frac{π}{2}"]) : ([-π / 2, 0.0, π / 2], ["", "", ""])

            trajs_ps = PlotSpec[]
            trajs = []
            foreach(u₀s[i], cs1) do (ν, γ), col
                u₀ = [ν, γ, 0.0, 0.0]
                sol = simulate(full_dynamics!, u₀, p, (0.0, 20.0))
                push!(trajs_ps, trajectory_phasespace(sol; color = col)...)
                push!(trajs, trajectory(sol, 0; color = col))
            end

            traj_axs = [
                (k, 1) => trajectory_axis(
                        [
                            trajs[k];
                        ]; xlabel = "", ylabel = "", xticksvisible = false, xticklabelsvisible = false, yticksvisible = false, yticklabelsvisible = false,
                    ) for k in eachindex(u₀s[i])
            ]

            push!(
                gls,
                (i, j) => S.GridLayout(
                    [
                        (1, 1) => phasespace_axis(
                            [
                                phasespace_regions(p; extend);
                                trajs_ps;
                                S.Scatter(νstar, γstar; marker = :xcross, color = :black, markersize = 12)
                            ]; extend, xlabel, ylabel, xticks, yticks, xticklabelspace = 5.0, yticklabelspace = 15.0
                        ),
                        (1, 1) => S.GridLayout(
                            S.Box(color = :white, width = Relative(0.4), height = Relative(0.4), halign = 0.97, valign = 0.97)
                        ),
                        (1, 1) => S.GridLayout(
                            diagram_axis(
                                [
                                    free_body_diagram(1.0, 0.3, p; pt = [0.55, 0.4], sf = 0.34, velocity = false, weight = false, gravity = false, labels = false, axes = false, angles = false, birdsize = 40, arrowheadsize = 10);
                                ]; backgroundcolor = :white
                            );
                            rowsizes = [Relative(0.4)], colsizes = [(Relative(0.4))], tellheight = false, tellwidth = false, halign = 0.97, valign = 0.97
                        ),
                        (1, 2) => S.GridLayout(traj_axs),
                        (1, 1) => S.Label(
                            text = axlabels1[i] * axlabels2[j], font = :bold, fontsize = 14pt, tellheight = false, halign = 0.0, tellwidth = false, padding = 0, valign = 1.17
                        ),
                    ],
                    colsizes = [Relative(0.6), Relative(0.4)]
                )
            )
        end
    end

    plot(f[1, 1], S.GridLayout(gls))

    rowgap!(f.layout.content[1].content, 18pt)
    colgap!(f.layout.content[1].content, 5)

    for i in eachindex(gls)
        rowgap!(f.layout.content[1].content.content[i].content.content[4].content, 5)
        colgap!(f.layout.content[1].content.content[i].content, 5)

        # add zoom inset for certain trajectories in physical space
        if i in [1, 2, 7, 8]
            idx = 7(i - 1) .+ [4, 5, 6] # indices of trajectory axis
            for ax in f.content[idx]
                zi = zoom_inset!(ax, Rect2f(-0.25, -2.5, 1.5, 3.0); inset_width = 0.35, inset_height = 0.7, halign = 0.95, valign = 0.7)
                hidedecorations!(zi.ax_inset)
            end
        end
    end

    return f
end
