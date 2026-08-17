using UnpoweredLandingsPaper, PointMass


## compare the accuracy of the direct collocation solutions to that with numerical integration
let
    # dc solution
    pa = TrigAerodynamics(1.5, 0.03, 2.0)
    u₀ = [1.5, 0.1, 1.0, 1.0]
    dc_sol = solve_ocp(dc_min_speed_unconstrained(u₀, pa); grid_size=200)

    # interpolate control and intergrate dynamics
    pc = InterpolatedControl(dc_sol.cntrlf)
    T = dc_sol.t[end]
    itp_sol = simulate(full_dynamics!, u₀, Parameters(pa, pc), (0.0, T); saveat=dc_sol.t)

    # percentage RMSE
    N = length(dc_sol.u)
    y1 = stack(itp_sol.u)
    y2 = stack(dc_sol.u)
    rmspe = sqrt.(mean(((y2 - y1)) .^ 2; dims=2))

    rmspe
end