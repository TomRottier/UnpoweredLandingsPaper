import CairoMakie.Makie.SpecApi as S

export trajectory, trajectory_phasespace, variable_timeseries, alpha_timeseries
export trajectory_axis, phasespace_axis, timeseries_axis, perching_axis
export phasespace_regions
export diagram_axis, stall_speed_diagram, term_speed_diagram, free_body_diagram
export force_vectors, velocity_vector
export add_subplot_labels

function equal_space(sol, n)
    return map(range(sol.u[1][3], sol.u[end][3], length = n)) do χ
        idx = findmin(x -> abs(x[3] - χ), sol.u)[2]
        return sol.t[idx]
    end
end


function timeseries_axis(plots = PlotSpec[], ylabel = ""; kwargs...)
    return S.Axis(;
        plots,
        xlabel = L"time $t / \hat{t}$", ylabel, limits = (0, nothing, nothing, nothing), kwargs...
    )
end

function variable_timeseries(sol::ODESolution, varidx; color = :gray, decorate = true, decmarkersize = 8, linestyle = :solid, kwargs...)
    plts = PlotSpec[]
    push!(plts, S.Lines(sol.t, sol[varidx, :]; color, linestyle, kwargs...))
    if decorate
        T = sol.t[end] - sol.t[begin]
        push!(
            plts,
            S.Scatter(Point2f(sol.t[1], sol[varidx, 1]); color, marker = :circle, markersize = decmarkersize, kwargs...),
            S.Scatter(Point2f(sol.t[end], sol[varidx, end]); color, marker = :diamond, markersize = decmarkersize, kwargs...)
        )
    end

    return plts
end

function alpha_timeseries(sol::ODESolution; color = :gray, decorate = true, decmarkersize = 8, linestyle = :solid, kwargs...)
    plts = PlotSpec[]
    f(t) = PointMass.alpha(sol(t), sol.prob.p, t)
    push!(plts, S.Lines(sol.t, f; color, linestyle, kwargs...))
    if decorate
        T = sol.t[end] - sol.t[begin]
        push!(
            plts,
            S.Scatter(Point2f(sol.t[1], f(sol.t[1])); color, marker = :circle, markersize = decmarkersize, kwargs...),
            S.Scatter(Point2f(sol.t[end], f(sol.t[end])); color, marker = :diamond, markersize = decmarkersize, kwargs...)
        )
    end

    return plts
end

function phasespace_axis(plots = PlotSpec[]; extend = sqrt(1 / 0.03), kwargs...)
    return S.Axis(;
        xlabel = L"speed $v / \hat{v}$", ylabel = L"flight path angle $\gamma$",
        limits = (0.0, extend, -π / 2, π / 2 + 0.05),
        xticks = ([0, extend], ["0", L"\nu_T^*"]), yticks = ([-π / 2, 0.0, π / 2], [L"-$\frac{π}{2}$", "0", L"\frac{π}{2}"]),
        xticklabelspace = 5.0, yticklabelspace = 15.0,
        kwargs...,
        plots
    )
end

function phasespace_regions(p; ν = true, γ = true, n = 200, extend = nothing, kwargs...)
    plts = PlotSpec[]
    yulim = π / 2 + 0.1
    yllim = -π / 2 - 0.1

    # acceleration nullcline
    if ν
        νs, γs = acceleration_nullcline(p; n)
        if !isnothing(extend)
            νs = [νs; extend]
            push!(γs, yllim)
        end
        push!(
            plts,
            S.Band(Point2.(νs, γs), Point2.(νs, yulim), color = cs3[1], alpha = 0.5),
        )
        push!(
            plts,
            S.Lines(νs, γs, color = cs3[2]),
            S.Band(Point2.(νs, yllim), Point2.(νs, γs), color = cs3[2], alpha = 0.5)
        )

    end

    # flight path nullcline
    if γ
        νs, γs = flightpath_nullcline(p; n)
        push!(
            plts,
            S.Lines([νs; reverse(νs)], [γs; -reverse(γs)], color = cs3[3]),
            S.Band(Point2.(νs[1:(end ÷ 2)], -γs[1:(end ÷ 2)]), Point2.(νs[1:(end ÷ 2)], γs[1:(end ÷ 2)]), color = cs3[3], alpha = 0.5)
        )
    end

    return plts
end

function trajectory_axis(plots = PlotSpec[]; kwargs...)
    return S.Axis(
        xlabel = L"range $x / \hat{l}$", ylabel = L"altitude $z / \hat{l}$",
        autolimitaspect = 1, xticklabelspace = 20.0, yticklabelspace = 20.0,
        xautolimitmargin = (0.1, 0.1);
        kwargs...,
        plots
    )
end


function trajectory(
        sol::ODESolution, nmarkers;
        color = :gray, decorate_start = true, decorate_end = true, evenspace = false, τs = nothing, forces = false,
        vel = false, weight = false, label = "", markersize = 30, decmarkersize = 8, sf = 1.0, kwargs...
    )
    plts = PlotSpec[]
    pts = [Point2(sol(t)[3], sol(t)[4]) for t in sol.t]
    l1 = S.Lines(pts; label, color, kwargs...)
    push!(plts, l1)

    if decorate_start
        push!(
            plts,
            S.Scatter(sol.u[1][3], sol.u[1][4]; marker = :circle, markersize = decmarkersize, color)
        )
        T = (sol.t[end] - sol.t[begin])
        T01 = sol.t[begin] + 0.1T
        push!(
            plts,
            S.Scatter(
                Point2(sol(T01)[3:4]), rotation = sol(T01)[2],
                marker = arrowhead, color = :black, markersize = 6
            ),
        )
    end

    if decorate_end
        push!(
            plts,
            S.Scatter(sol.u[end][3], sol.u[end][4]; marker = :diamond, markersize = decmarkersize, color)
        )
    end

    if nmarkers !== 0
        if evenspace && isnothing(τs)
            τs = equal_space(sol, nmarkers)
        elseif isnothing(τs)
            τs = range(sol.t[1], sol.t[end], nmarkers)
        end
        pts = map(τ -> Point2(sol(τ; idxs = 3:4)...), τs)
        rots = map(τ -> sol(τ, idxs = 2) + alpha(sol(τ), sol.prob.p, τ), τs)

        if forces
            foreach(pts, τs) do pt, τ
                push!(
                    plts,
                    force_vectors(pt, sol(τ, idxs = 1:2)..., sol.prob.p; sf, weight, drawbird = false)...
                )
            end
        end

        if vel
            foreach(pts, τs) do pt, τ
                push!(plts, velocity_vector(pt, [sf, 1] .* sol(τ, idxs = 1:2)...)...)
            end
        end

        push!(
            plts,
            S.Scatter(pts; rotation = rots, marker = bird, markersize, color, label, filter(!=(:linestyle) ∘ first, kwargs)...)
        )
    end

    return plts
end

function trajectory(
        sol::DirectCollocationSolution, nmarkers;
        color = :gray, decorate_start = true, decorate_end = true, evenspace = false, τs = nothing, decmarkersize = 8, label = "", markersize = 30, kwargs...
    )
    plts = PlotSpec[]
    pts = [Point2(u[3:4]) for u in sol.u]
    l1 = S.Lines(pts; label, color, kwargs...)
    push!(plts, l1)

    if decorate_start
        push!(
            plts,
            S.Scatter(sol.u[1][3], sol.u[1][4]; marker = :circle, markersize = decmarkersize, color)
        )
        T = (sol.t[end] - sol.t[begin])
        T01 = sol.t[begin] + 0.1T
        push!(
            plts,
            S.Scatter(
                Point2(sol(T01)[3:4]), rotation = sol(T01)[2],
                marker = arrowhead, color = :black, markersize = 6
            ),
        )
    end

    if decorate_end
        push!(
            plts,
            S.Scatter(sol.u[end][3], sol.u[end][4]; marker = :diamond, markersize = decmarkersize, color)
        )
    end

    if nmarkers !== 0
        if evenspace && isnothing(τs)
            τs = equal_space(sol, nmarkers)
        elseif isnothing(τs)
            τs = range(sol.t[1], sol.t[end], nmarkers)
        end
        pts = map(τ -> Point2(sol(τ)[3:4]), τs)
        rots = map(τ -> sol(τ)[2] + sol.cntrlf(τ)[1], τs)

        push!(
            plts,
            S.Scatter(pts; rotation = rots, marker = bird, markersize, color, label, filter(!=(:linestyle) ∘ first, kwargs)...)
        )
    end

    return plts
end

function force_vectors(
        pt, ν, γ, p::Parameters{<:Aerodynamics, <:Control};
        sf = 1.0, weight = false, arrowheadsize = 15, decmarkersize = 8, kwargs...
    )
    plts = PlotSpec[]

    α = alpha([ν, γ], p, nothing)
    cl = CL(α, p)
    cd = CD(α, p)
    W = sf * -1.0
    L = sf * cl * ν^2
    D = sf * cd * ν^2

    push!(
        plts,
        S.Arrows2D(pt, Point2(-L * sin(γ), L * cos(γ)); color = lift_colour, tip = arrowhead_force, tipwidth = arrowheadsize, tiplength = arrowheadsize, kwargs...),
        S.Arrows2D(pt, Point2(-D * cos(γ), -D * sin(γ)); color = drag_colour, tip = arrowhead_force, tipwidth = arrowheadsize, tiplength = arrowheadsize, kwargs...)
    )

    if weight
        push!(
            plts,
            S.Arrows2D(pt, Point2(0.0, W); color = weight_colour, tip = arrowhead_force, tipwidth = arrowheadsize, tiplength = arrowheadsize, kwargs...)
        )
    end

    return plts
end

function velocity_vector(pt, ν, γ; sf = 1.0, color = :black, arrowheadsize = 8, kwargs...)
    ν *= sf
    νx = ν * cos(γ)
    νz = ν * sin(γ)
    pt2 = pt .+ [νx, νz]
    return [
        S.Lines(Point2f[pt, pt2]; color, linewidth = 2, kwargs...),
        S.Scatter(Point2f(pt2); rotation = γ, color, marker = Polygon(arrowhead_vel), markersize = arrowheadsize, strokewidth = 0, kwargs...),
    ]
    # S.Arrows2D(pt, Point2(ν * cos(γ), ν * sin(γ)); color, tip=arrowhead_vel, tipwidth=arrowheadsize, tiplength=arrowheadsize, kwargs...)
end


function trajectory_phasespace(
        sol::ODESolution;
        color = :gray, decorate_start = true, decorate_end = true, decmarkersize = 8, kwargs...
    )
    plts = PlotSpec[]
    pts = [Point2(sol(t)[1], sol(t)[2]) for t in sol.t]
    l1 = S.Lines(pts; color, kwargs...)
    push!(plts, l1)

    T0, T = sol.t[begin], (sol.t[end] - sol.t[begin])
    if decorate_start
        T01 = T0 + 0.1T
        push!(
            plts,
            S.Scatter(Point2(sol(T0, idxs = 1:2)); marker = :circle, color, markersize = decmarkersize),
            S.Scatter(Point2(sol(T01)[1:2]), rotation = atan(sol(T01, Val{1})[[2, 1]]...), marker = arrowhead, color = :black, markersize = 6)
        )
    end
    if decorate_end
        push!(
            plts,
            S.Scatter(Point2(sol(T + T0, idxs = 1:2)); marker = :diamond, color, markersize = decmarkersize)
        )
    end

    return plts
end


function diagram_axis(plots = PlotSpec[]; kwargs...)
    return S.Axis(
        plots = plots,
        autolimitaspect = 1,
        xlabelvisible = false,
        ylabelvisible = false,
        xticklabelsvisible = false,
        yticklabelsvisible = false,
        xticksvisible = false,
        yticksvisible = false;
        kwargs...
    )
end

function stall_speed_diagram(ν, γ, p; x = 0.5, y = 0.5, sf = 1.0, narrows = 9, textx = 0.7, texty = 0.15, arrow_alpha = 1.0, forces = true, fontsize = 12pt, markersize = 70, markercolor = :gray, fbd_arrows = true)
    plts = PlotSpec[]
    space = :relative
    α = PointMass.alpha([ν, γ], p, nothing)
    cl = CL(α, p)

    # stall speed
    γs = range(-π / 2, π / 2, length = narrows)
    νxs = [sf * sqrt(cos(γ)^3 / cl) + x for γ in range(-π / 2, π / 2, 101)]
    νzs = [
        [sf * sqrt(cos(γ) * sin(γ)^2 / cl) + y for γ in range(-π / 2, 0.0, 50)];
        [sf * -sqrt(cos(γ) * sin(γ)^2 / cl) + y for γ in range(0.0, π / 2, 51)]
    ]
    push!(plts, S.Lines(νxs, νzs; color = cs3[3], space))

    foreach(γs) do γ
        _ν = sqrt(cos(γ) / cl)
        νx = _ν * cos(γ) * sf
        νz = _ν * sin(γ) * sf
        L = cl * _ν^2 * sf
        Lx = -L * sin(γ)
        Lz = L * cos(γ)
        push!(
            plts,
            S.Lines(Point2f[(x, y), (νx + x, νz + y)]; color = :gray, linewidth = 2, alpha = arrow_alpha, space),
            S.Scatter(Point2(νx + x, νz + y); marker = Polygon(arrowhead_vel), markersize = 7, color = :gray, strokewidth = 0, rotation = γ, alpha = arrow_alpha, space),
        )
        if forces
            push!(
                plts,
                S.Arrows2D(Point2(νx + x, νz + y), Point2(Lx, Lz); color = lift_colour, alpha = arrow_alpha, space)
            )
        end
    end

    push!(
        plts,
        free_body_diagram(
            ν, γ, p;
            pt = [x, y], sf, axes = false, angles = false, freestream = false, labels = false, gravity = false, velocity = fbd_arrows, forces = fbd_arrows, weight = fbd_arrows, birdsize = markersize, markercolor
        )...,
        S.Text(textx, texty; text = L"\nu_{\text{stall}} = \sqrt{\frac{\cos \gamma}{C_L}}", space, align = (:center, :center), color = cs3[3], fontsize)
    )

    return plts
end

function term_speed_diagram(ν, γ, p; x = 0.5, y = 0.5, sf = 1.0, narrows = 9, textx = 0.65, texty = 0.75, arrow_alpha = 1.0, forces = true, fontsize = 12pt, markersize = 70, markercolor = :gray, fbd_arrows = true)
    plts = PlotSpec[]
    space = :relative
    α = PointMass.alpha([ν, γ], p, nothing)
    cd = CD(α, p)

    # terminal speed
    γs = range(-π / 2, 0.0, narrows)
    νxs = [sf * sqrt(-sin(γ) * cos(γ)^2 / cd) + x for γ in range(0.0, -π / 2, 101)]
    νzs = [sf * -sqrt(-sin(γ)^3 / cd) + y for γ in range(0.0, -π / 2, 101)]
    push!(plts, S.Lines(νxs, νzs; color = cs3[2], space))

    foreach(γs) do γ
        _ν = sqrt(-sin(γ) / cd)
        νx = _ν * cos(γ) * sf
        νz = _ν * sin(γ) * sf
        D = cd * _ν^2 * sf
        Dx = -D * cos(γ)
        Dz = -D * sin(γ)
        push!(
            plts,
            S.Lines(Point2f[(x, y), Point2(νx + x, νz + y)]; color = :gray, linewidth = 2, alpha = arrow_alpha, space),
            S.Scatter(Point2(νx + x, νz + y); marker = Polygon(arrowhead_vel), markersize = 7, color = :gray, strokewidth = 0, rotation = γ, alpha = arrow_alpha, space),
        )
        if forces
            push!(
                plts,
                S.Arrows2D(Point2(νx + x, νz + y), Point2(Dx, Dz); color = drag_colour, alpha = arrow_alpha, space)
            )
        end
    end

    # draw bird
    return push!(
        plts,
        free_body_diagram(
            ν, γ, p;
            pt = [x, y], sf, axes = false, angles = false, freestream = false, labels = false, gravity = false, velocity = fbd_arrows, forces = fbd_arrows, weight = fbd_arrows, birdsize = markersize, markercolor
        )...,
        S.Text(textx, texty; text = L"\nu_{\text{term}} = \sqrt{\frac{\minus \sin \gamma}{C_D}}", space, align = (:center, :center), color = cs1[2], fontsize)
    )


end

#### MADE FOR FIGURE SIZE (350,350) - USING VERY DIFFERENT SIZES WILL MESS THINGS UP ####
function free_body_diagram(
        ν, γ, p, traj = nothing; pt = [0.5, 0.5], gpt = [0.8, 0.3], upt = [0.95, 0.8],
        sf = 1.0, axes = true, angles = true, velocity = true, forces = true, weight = true, gravity = true, labels = true, freestream = true, birdsize = 90, arrowheadsize = 15, markercolor = :gray
    )
    plts = PlotSpec[]

    space = :relative
    α = PointMass.alpha([ν, γ], p, nothing)
    θ = γ + α
    _ν = ν * sf
    L = CL(α, p) * ν^2 * sf
    D = CD(α, p) * ν^2 * sf


    if !isnothing(traj)
        # rescale trajectory to begin at `pt` and in relative units [0,1]
        x0, x1 = traj[1][1], traj[end][1]
        traj_scaled = map(p -> Point2((p) / (x1) + pt), traj)
        push!(plts, S.Lines(traj_scaled; color = to_color((0.7, 0.7, 0.7)), space))
    end

    # angles
    if angles
        push!(
            plts,
            S.Lines(Point2f[pt, pt + [_ν, 0.0]]; linestyle = (:dash, :dense), color = :black, linewidth = 1, space),
            S.Lines(Point2f[pt, pt + _ν * [cos(θ), sin(θ)]]; linestyle = (:dash, :dense), color = :black, linewidth = 1, space),
            S.Lines(
                [Point2(pt + 0.5 * _ν * [cos(ψ), sin(ψ)]) for ψ in range(0.0, γ, 50)];
                color = :black, linewidth = 1, space
            ),
            S.Scatter(
                Point2(pt + 0.5 * _ν * [cos(γ), sin(γ)]);
                marker = Polygon(arrowhead_vel), rotation = γ + π / 2, markersize = 5, color = :black, space
            ),
            S.Lines(
                [Point2(pt + 0.75 * _ν * [cos(ψ), sin(ψ)]) for ψ in range(γ, θ, 50)];
                color = :black, linewidth = 1, space
            ),
            S.Scatter(
                Point2(pt + 0.75 * _ν * [cos(θ), sin(θ)]);
                marker = Polygon(arrowhead_vel), rotation = θ + π / 2, markersize = 5, color = :black, space
            ),
        )
    end

    # velocity
    if velocity
        push!(
            plts,
            # bird + vectors
            velocity_vector(Point2(pt), ν, γ; sf, arrowheadsize = 0.5arrowheadsize, space)...,
        )
    end

    # forces
    if forces
        push!(
            plts,
            force_vectors(Point2(pt), ν, γ, p; sf, weight, arrowheadsize, space)...,
        )
    end

    # gravity
    if gravity
        push!(
            plts,
            S.Lines(Point2f[(gpt - [0.0, 0.1]), (gpt + [0.0, 0.1])]; color = :black, linewidth = 2, space),
            S.Scatter(
                Point2f[gpt, gpt];
                marker = Polygon(arrowhead_vel), rotation = -π / 2, color = :black, markersize = 7,
                marker_offset = [(0, -5), (0, 5)], space
            ),
        )
    end

    # bird
    push!(
        plts,
        S.Scatter(Point2(pt); marker = bird, color = markercolor, markersize = birdsize, rotation = θ, space)
    )

    # freestream
    if freestream
        push!(
            plts,
            velocity_vector(Point2(upt), ν, γ + π; sf, arrowheadsize = 0.5arrowheadsize, space)...,
        )
    end

    # labels
    if labels
        push!(
            plts,
            S.Text(
                Point2(pt + 0.65 * _ν * [cos(γ / 2), sin(γ / 2)]);
                text = L"\gamma", align = (:center, :center), space
            ),
            S.Text(
                Point2(pt + 0.9 * _ν * [cos(γ + α / 2), sin(γ + α / 2)]);
                text = L"\alpha", align = (:center, :center), space
            ),
            S.Text(
                Point2(pt + L * [-sin(γ), cos(γ)]);
                text = L"L", align = (:center, :top), offset = (-10, -10), color = lift_colour, space
            ),
            S.Text(
                Point2(pt - D * [cos(γ), sin(γ)]);
                text = L"D", align = (:center, :center), offset = (0, 15), color = drag_colour, space
            ),
            S.Text(
                Point2(pt + [0, -sf]);
                text = L"m g", align = (:center, :bottom), offset = (15, 0), space
            ),
            S.Text(
                Point2(pt + _ν * [cos(γ), sin(γ)]);
                text = L"v", align = (:center, :center), offset = (0, 8), space
            ),
            S.Text(
                Point2(gpt);
                text = L"g", align = (:center, :center), offset = (15, 0), space
            ),
            S.Text(
                Point2(upt);
                text = L"U_{\infty}", align = (:center, :center), offset = (-25, 5), space
            )
        )
    end

    # axes
    if axes
        axpt = [0.15, 0.2]
        axlen = 0.19
        lo = 9 # label offset
        push!(
            plts,
            S.Arrows2D(
                Point2(axpt), Point2(axlen, 0.0);
                tip = arrowhead_axis, tipwidth = 12, tiplength = 12, space
            ),
            S.Arrows2D(
                Point2(axpt), Point2(0.0, -axlen);
                tip = arrowhead_axis, tipwidth = 12, tiplength = 12, space
            ),
            S.Arrows2D(
                Point2(axpt), Point2(axlen * cos(γ), axlen * sin(γ));
                tip = arrowhead_axis, tipwidth = 12, tiplength = 12, space
            ),
            S.Arrows2D(
                Point2(axpt), Point2(-axlen * sin(γ), axlen * cos(γ));
                tip = arrowhead_axis, tipwidth = 12, tiplength = 12, space
            ),
            S.Scatter(Point2(axpt); marker = :circle, color = :white, markersize = 15, space),
            S.Scatter(Point2(axpt); marker = :circle, color = :black, markersize = 3, space),
            S.Lines(
                [Point2(axpt + 0.5 * axlen * [cos(ψ), sin(ψ)]) for ψ in range(0.0, γ, 50)];
                color = :black, linewidth = 1, space
            ),
            S.Scatter(
                Point2(axpt + 0.5 * axlen * [cos(γ), sin(γ)]);
                marker = Polygon(arrowhead_vel), rotation = γ + π / 2, markersize = 5, color = :black, space
            ),
            S.Text(
                Point2(axpt .+ [axlen, 0.0]);
                text = L"\hat{e}_x", align = (:center, :center), offset = (lo, 0), space
            ),
            S.Text(
                Point2(axpt .+ [0.0, -axlen]);
                text = L"\hat{e}_z", align = (:center, :center), offset = (15, 8), space
            ),
            S.Text(
                Point2(axpt .+ axlen * [cos(γ), sin(γ)]);
                text = L"\hat{e}_v", align = (:center, :center), offset = lo .* ((cos(γ)), sin(γ)), space
            ),
            S.Text(
                Point2(axpt .+ 1.2axlen * [-sin(γ), cos(γ)]);
                text = L"\hat{e}_\gamma", align = (:center, :center), offset = (10, 0), space
            ),
            S.Text(
                Point2(axpt .+ 0.65 * axlen .* [cos(γ / 2), sin(γ / 2)]);
                text = L"\gamma", align = (:center, :center), space
            ),
        )
    end

    return plts
end


function perching_axis(plots = PlotSpec[]; skip = 1, perch_width = 0.1, data = true, kwargs...)
    ph = experimental_data.perch_height
    pw = perch_width * ph
    ps = experimental_data.perch_spacing

    plts = [
        S.HLines(0.0, color = :black, linewidth = 1),
        S.Lines(
            [-ps - pw, -ps, -ps + pw, NaN, -pw, 0.0, pw], [0.0, ph, 0.0, NaN, 0.0, ph, 0.0];
            color = :black, linewidth = 2
        ),
    ]

    if data
        χ = experimental_data.χ[11:skip:end] # defined start based on Kleinheerenbrink et al.
        ζ = experimental_data.ζ[11:skip:end]
        push!(
            plts,
            S.Scatter(χ, ζ; color = :black, marker = :circle, markersize = 6)
        )
    end

    push!(plts, plots...)

    return trajectory_axis(plts; kwargs...)
end
