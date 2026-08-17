export simulate
export acceleration_nullcline, flightpath_nullcline

isoutofdomain(u, p, t) = false

"""
    simulate(f, u0, p, tspan=(0.0, 5.0), tspan=(0.0, 5.0); kwargs...)

Numerically integrate f starting from u0.

`kwargs` get passed on to the call to `solve`

default tolerances of `abstol=1e-6` and `reltol=1e-6` but can be overidden with in `kwargs`
"""
function simulate(f, u0, p, tspan=(0.0, 5.0); kwargs...)
    # if simulating backwards in time, repeat simulation forwards to make solution work with other functions (expecting time increasing)
    if tspan[1] > tspan[2]
        prob = ODEProblem(f, u0, tspan, p; save_everystep=false, abstol=1.0e-6, reltol=1.0e-6)
        sol = solve(prob, Tsit5(); isoutofdomain, kwargs...)
        newprob = remake(prob, u0=sol.u[end], tspan=(0.0, -sol.t[end]), saveat=0.01, save_end=false)
        kwargs = filter(!=(:callback) ∘ first, kwargs) # remove callback on second call
        return solve(newprob, Tsit5(); isoutofdomain, kwargs...)
    end

    prob = ODEProblem(f, u0, tspan, p; saveat=0.01, abstol=1.0e-6, reltol=1.0e-6)
    return solve(prob, Tsit5(); isoutofdomain, kwargs...)
end

"""
    acceleration_nullcline(cd, p; n=100)

Returns the nullcline for ν `(νs, γs)` for a constant `cd`
"""
function acceleration_nullcline(p::Parameters{<:Aerodynamics,ConstantAlpha}; n=100)
    cd = CD(alpha(nothing, p, nothing), p)
    νs = range(0.0, min(sqrt(1 / cd), 5.77), length=n)
    γs = [asin(max(-cd * ν^2, -1.0)) for ν in νs]
    return (νs, γs)
end

"""
    flightpath_nullcline(cl, p; n=100)

Returns the nullcline for γ `(νs, γs)` for a constant `cl`. Symmetrical about γ = 0.
"""
function flightpath_nullcline(p::Parameters{<:Aerodynamics,ConstantAlpha}; n=100)
    cl = CL(alpha(nothing, p, nothing), p)
    cd = CD(alpha(nothing, p, nothing), p)
    νs = range(0.0, min(sqrt(1 / cl), 5.8), length=n ÷ 2)
    γs = [acos(min(cl * ν^2, 1.0)) for ν in νs]
    return ([νs; reverse(νs)], [-γs; reverse(γs)])
end
