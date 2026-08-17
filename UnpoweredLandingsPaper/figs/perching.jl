export perching

function compare_post_stall_time_to_dist(pa, CTmax, stall_angle)
    init = (time = [0.0, 3.0], control = [[0.0, CTmax], [π / 2, 0.0]])
    sol0 = solve_ocp(dc_weighted_min_speed_poststall_time(pa, 1.0; CTmax, stall_angle); method = :exa, display = false, init)
    sol1 = solve_ocp(dc_weighted_min_speed_poststall_time(pa, 0.0; CTmax, stall_angle); method = :exa, display = false, init)
    sol2 = solve_ocp(dc_weighted_min_speed_poststall_dist(pa, 0.0; CTmax, stall_angle); method = :exa, display = false, init)

    sols = [sol0, sol1, sol2]
    cols = [1.0, 0.0, :purple]

    ax0 = perching_axis(
        map(sols, cols) do sol, col
            return S.Lines([Point2f(u[3:4]) for u in sol.u]; color = col, colorrange = (0.0, 1.0))
        end;
        data = false, yticks = [0.0, 0.2, 0.4], xticklabelspace = 10.0
    )

    ax2a = timeseries_axis(
        map(sols, cols) do sol, col
            νs = [sol(t)[1] for t in sol.t]
            S.Lines(sol.t, νs; color = col, colorrange = (0.0, 1.0))
        end,
        L"v / \hat{v}"; yticks = [0.4, 0.8, 1.2]
    )
    ax2b = timeseries_axis(
        mapreduce(vcat, sols, cols) do sol, col
            αs = [sol.cntrlf(t)[1] for t in sol.t]
            [
                S.HLines(stall_angle; color = :gray, linestyle = (:dash, :dense));
                S.Lines(sol.t, αs; color = col, colorrange = (0.0, 1.0))
            ]
        end,
        L"\alpha"; yticks = ([0, π / 4, π / 2], ["0", L"\frac{π}{4}", L"\frac{π}{2}"])
    )

    leg = S.Legend(
        [LineElement(color = col, colorrange = (0.0, 1.0), linewidth = 3) for col in cols],
        ["speed", "time", "distance"];
        orientation = :horizontal
    )

    return S.GridLayout(
        (0, 1) => leg,
        (1, 1) => ax0,
        (2, 1) => S.GridLayout([ax2a ax2b])
    )

end

function cost_weighting(pa, CTmax, stall_angle; n = 11)
    λs = range(0.0, 1.0, length = n)
    init = nothing
    sols = Folds.map(λs) do λ
        solve_ocp(dc_weighted_min_speed_poststall_time(pa, λ; CTmax, stall_angle); method = :exa, display = false, init)
    end

    return S.GridLayout(
        perching_axis(
            mapreduce(vcat, sols, λs) do sol, λ
                alpha = λ ∈ (0.0, 1.0) ? 1.0 : 0.3
                bird_traj = λ ∈ (0.0, 1.0) ? trajectory(sol, 4; evenspace = true, colorrange = extrema(λs), color = λ, alpha) : []
                [
                    S.Lines([Point2f(u[3:4]) for u in sol.u]; colorrange = extrema(λs), color = λ, alpha);
                    bird_traj
                ]

            end;
            skip = 3, yticklabelspace = 25.0, title = L"C_{T\text{max}} = %$CTmax"
        )
    )
end

function touchdown_speeds(pa, _CTmax, stall_angle; n = 11)
    CTmaxs = _CTmax isa AbstractArray ? _CTmax : [_CTmax]
    λs = range(1.0, 0.0, length = n)
    νtds = map(CTmaxs) do CTmax
        sol0 = solve_ocp(dc_weighted_min_speed_poststall_time(pa, 0.0; CTmax, stall_angle); method = :exa, display = false)
        sol1 = solve_ocp(dc_weighted_min_speed_poststall_time(pa, 1.0; CTmax, stall_angle); method = :exa, display = false)
        νmin, νmax = objective(sol1.sol), sol0.u[end][1]
        Tmin, Tmax = objective(sol0.sol), quadgk(t -> post_stall_x(1.0, sol1.cntrlf(t)[1], stall_angle), sol1.t[1], sol1.t[end])[1]
        init = sol1.sol
        return Folds.map(λs) do λ
            sol = solve_ocp(dc_weighted_min_speed_poststall_time(pa, λ; CTmax, stall_angle, νmin, νmax, Tmin, Tmax); method = :exa, display = false, init)
            init = sol.sol
            return sol.u[end][1]
        end
    end

    return S.GridLayout(
        S.Axis(
            plots = [
                [S.Lines(λs, νtd; color = λs) for νtd in νtds];
                S.Text(Point2(1.0, νtds[1][1]); text = L"C_{T\text{max}} = %$(CTmaxs[1])", align = (:right, :center), offset = (-5, 13));
                S.Text(Point2(0.6, νtds[end][1]); text = L"C_{T\text{max}} = %$(CTmaxs[end])", align = (:right, :center), offset = (0, -2))
            ],
            xlabel = L"\lambda", ylabel = L"touchdown speed $v / \hat{v}$"
        )
    )
end


function perching()
    f = Figure(size = (full_page_fig_width, 0.8full_page_fig_height), figure_padding = (5, 10, 5, 5))
    pa = TrigAerodynamics(1.5, 0.03, 2.0)
    stall_angle = deg2rad(20)

    CTmaxs = [0.8, 1.0, 1.2, 1.4]
    gls = map(CTmaxs) do CTmax
        cost_weighting(pa, CTmax, stall_angle)
    end
    gls[1].content[1].second.xlabelvisible = false
    gls[1].content[1].second.xticklabelsvisible = false
    gls[2].content[1].second.xlabelvisible = false
    gls[2].content[1].second.xticklabelsvisible = false
    gls[2].content[1].second.ylabelvisible = false
    gls[2].content[1].second.yticklabelsvisible = false
    gls[4].content[1].second.ylabelvisible = false
    gls[4].content[1].second.yticklabelsvisible = false


    cb1a = S.Colorbar(; limits = (0.0, 1.0), label = "cost weighting, λ", vertical = false, ticks = ([0.0], ["post-stall time"]), ticklabelalign = (:left, :bottom))
    cb1b = S.Colorbar(; limits = (0.0, 1.0), label = "cost weighting, λ", vertical = false, ticks = ([1.0], ["touchdown speed"]), ticklabelalign = (:right, :bottom))
    cb2 = S.Colorbar(; limits = (0.0, 1.0), vertical = false, ticks = [0.0, 1.0], flipaxis = false)


    gl1 = S.GridLayout(
        [
            (1, 1, TopLeft()) => S.Label(; text = "C", padding = (0, 0, 10, 0), font = :bold, fontsize = 14pt),
            (1, 1) => gls[1],
            (1, 2, TopLeft()) => S.Label(; text = "D", padding = (0, 0, 10, 0), font = :bold, fontsize = 14pt),
            (1, 2) => gls[2],
            (2, 1, TopLeft()) => S.Label(; text = "E", padding = (0, 0, 10, 0), font = :bold, fontsize = 14pt),
            (2, 1) => gls[3],
            (2, 2, TopLeft()) => S.Label(; text = "F", padding = (0, 0, 10, 0), font = :bold, fontsize = 14pt),
            (2, 2) => gls[4],
        ]
    )

    gl2 = touchdown_speeds(pa, [CTmaxs;], stall_angle; n = 20)

    gl3 = compare_post_stall_time_to_dist(pa, 0.8, stall_angle)

    gl = S.GridLayout(
        [
            (0, :) => cb1a,
            (0, :) => cb1b,
            (0, :) => cb2,
            (1, 1) => S.GridLayout([(1, 1) => gl3, (1, 1, TopLeft()) => S.Label(; text = "A", font = :bold, fontsize = 14pt)]),
            (1, 2) => S.GridLayout([(1, 1) => gl2, (1, 1, TopLeft()) => S.Label(; text = "B", font = :bold, fontsize = 14pt)]),
            (2, :) => gl1,
        ]
    )

    plot(f[1, 1], gl)
    linkaxes!(f.content[(end - 6):2:end]...)
    rowgap!(f.layout.content[1].content, 0)

    return f
end
