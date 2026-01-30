export DirectCollocationSolution, solve_ocp
export dc_min_speed_unconstrained
export dc_min_speed_ground, dc_min_speed_constrain_td_angle
export dc_weighted_min_speed_poststall_dist, dc_weighted_min_speed_poststall_time

struct DirectCollocationSolution{T1,T2,F1,F2}
    sol::T1
    t::Vector{T2}
    u::Vector{Vector{T2}}
    statef::F1
    cntrlf::F2

    function DirectCollocationSolution(sol::T1) where {T1}
        ts = time_grid(sol)
        u = state(sol)
        c = control(sol)
        us = [u(t) for t in ts]
        return new{T1,eltype(ts),typeof(u),typeof(c)}(sol, ts, us, u, c)
    end
end
(sol::DirectCollocationSolution)(t) = sol.statef(t)


"""
    solve_ocp(ocp; inits=[nothing], display=true, disc_method=:midpoint, grid_size=100; kwargs...)

Solves the optimal control problem `ocp`. `inits` is a vector of named tuples specifying different initial guesses.

Procedure:
    - solve `ocp` for each init in `inits`, solved with grid_size = 25, disc_method = :trapeze
    - take solution with lowest objective
    - use this solution as initial guess for solve with grid_size=`grid_size` and disc_method = `disc_methods`

`init` allows first optimisation to be initialised.

"""
function solve_ocp(ocp; method=:adnlp, solver=:ipopt, disc_method=:midpoint, init=nothing, grid_size=100, display=true, kwargs...)
    retry_count = 0
    while retry_count < 20
        try
            sol = solve(ocp, :direct, method, solver; disc_method, grid_size=25, display=false, init)
            sol = solve(ocp, :direct, method, solver; disc_method, grid_size, display, init=sol)
            return DirectCollocationSolution(sol)

        catch e
            if e isa DomainError
                retry_count += 1
                @warn "DomainError, retrying (count $retry_count)"
            else
                rethrow(e)
            end
        end

    end

    error("solving optimal control problem resulted in too many reattempts due to DomainError's ($retry_count)")
end


"""
    lb(value, tol)

Sets lower bound of `value` for a given tolerance `tol` e.g. lb(1.0, 0.05) = 0.95. Correctly accounts for negative values.

"""
lb(val, tol) = min((1 - tol) * val, (1 + tol) * val)

"""
    ub(value, tol)

Sets upper bound of `value` for a given tolerance `tol` e.g. lb(1.0, 0.05) = 1.05. Correctly accounts for negative values.
    
"""
ub(val, tol) = max((1 - tol) * val, (1 + tol) * val)


"""
    dc_min_speed_unconstrained(u₀, pa)

Create optimal control problem which minimises touchdown speed from given initial conditions `u₀` and aerodynamic model `pa`. Otherwise completely unconstrained.

"""
function dc_min_speed_unconstrained(u₀, pa)
    return @def begin
        T ∈ R, variable
        t ∈ [0, T], time
        x = (ν, γ, χ, ζ) ∈ R⁴, state
        α ∈ R, control
        ν̇ = -CD(α(t), pa) * ν(t)^2 - sin(γ(t))
        γ̇ = CL(α(t), pa) * ν(t) - cos(γ(t)) / ν(t)
        ẋ(t) == [ν̇, γ̇, ν(t) * cos(γ(t)), ν(t) * sin(γ(t))]
        ν(t) ≥ 0.01
        -π / 2 ≤ γ(t) ≤ π / 2
        0.0 ≤ α(t) ≤ π / 2
        ν(0) == u₀[1]
        γ(0) == u₀[2]
        χ(0) == u₀[3]
        ζ(0) == u₀[4]
        ν̇ ≤ 0.0
        ν(T) → min
    end
end


"""
    dc_min_speed_constrain_td_angle(u₀, pa, γtd; Tmax=5.0)

Create optimal control problem which minimises touchdown speed from given initial conditions `u₀` and aerodynamic model `pa`. The touchdown speed is constrained to be `γtd`.

The vertical coordinate is constrained to finish at 0 for plotting purposes, ζ(0) in u₀ is ignored, -π/2 ≤ γ(t) ≤ 0.0, `Tmax` constrains the final time.

"""
function dc_min_speed_constrain_td_angle(u₀, pa, γtd; Tmin=0.0, Tmax=5.0, ν̇max=0.0)
    return @def begin
        T ∈ R, variable
        t ∈ [0, T], time
        x = (ν, γ, χ, ζ) ∈ R⁴, state
        α ∈ R, control
        ν̇ = -CD(α(t), pa) * ν(t)^2 - sin(γ(t))
        γ̇ = CL(α(t), pa) * ν(t) - cos(γ(t)) / ν(t)
        Tmin ≤ T ≤ Tmax
        ∂(ν)(t) == ν̇
        ∂(γ)(t) == γ̇
        ∂(χ)(t) == ν(t) * cos(γ(t))
        ∂(ζ)(t) == ν(t) * sin(γ(t))
        ν(t) ≥ 0.01
        -π / 2 ≤ γ(t) ≤ 0.0
        0.0 ≤ α(t) ≤ π / 2
        ν(0) == u₀[1]
        γ(0) == u₀[2]
        χ(0) == u₀[3]
        ζ(T) == 0.0
        γ(T) == γtd
        ν̇ ≤ ν̇max
        ν(T) → min
    end
end

# T ∈ R, variable
# t ∈ [0, T], time
# x = (ν, γ, χ, ζ, α) ∈ R⁵, state
# α̇ ∈ R, control
# ν̇ = -CD(α(t), pa) * ν(t)^2 - sin(γ(t))
# γ̇ = CL(α(t), pa) * ν(t) - cos(γ(t)) / ν(t)
# Tmin ≤ T ≤ Tmax
# ∂(ν)(t) == ν̇
# ∂(γ)(t) == γ̇
# ∂(χ)(t) == ν(t) * cos(γ(t))
# ∂(ζ)(t) == ν(t) * sin(γ(t))
# ∂(α)(t) == α̇(t)
# ν(t) ≥ 0.01
# -π / 2 ≤ γ(t) ≤ 0.0
# 0.0 ≤ α(t) ≤ π / 2
# -2π ≤ α̇(t) ≤ 2π
# ν(0) == u₀[1]
# γ(0) == u₀[2]
# χ(0) == u₀[3]
# ζ(T) == 0.0
# ν̇ ≤ ν̇max
# (ν(T)) → min


function dc_min_speed_ground(u₀, pa; Tmin=0.0, Tmax=5.0, ν̇max=0.0)
    return @def begin
        T ∈ R, variable
        t ∈ [0, T], time
        x = (ν, γ, χ, ζ) ∈ R⁴, state
        α ∈ R, control
        ν̇ = -CD(α(t), pa) * ν(t)^2 - sin(γ(t))
        γ̇ = CL(α(t), pa) * ν(t) - cos(γ(t)) / ν(t)
        Tmin ≤ T ≤ Tmax
        ∂(ν)(t) == ν̇
        ∂(γ)(t) == γ̇
        ∂(χ)(t) == ν(t) * cos(γ(t))
        ∂(ζ)(t) == ν(t) * sin(γ(t))
        ν(t) ≥ 0.01
        -π / 2 ≤ γ(t) ≤ 0.0
        0.0 ≤ α(t) ≤ π / 2
        ν(0) == u₀[1]
        γ(0) == u₀[2]
        χ(0) == u₀[3]
        ζ(T) == 0.0
        ν̇ ≤ ν̇max
        ν(T) → min
    end
end


post_stall_x(x, α, αstall) = x * 0.5 * (tanh((α - αstall) / deg2rad(1)) + 1)

"""
dc_weighted_min_speed_poststall_dist(p, λ; CTmin=0.0, CTmax=1.0, pos_bounds=0.05, stall_angle=deg2rad(25))

Create optimal control problem which minimises weighted cost between touchdown speed and post-stall *distance*. Weighted by parameter `λ`

`pos_bounds` chooses the tolerance to the measured touchdown position from the experimental data.
"""
function dc_weighted_min_speed_poststall_dist(p, λ; CTmin=0.0, CTmax=1.0, pos_bounds=0.05, stall_angle=deg2rad(25))
    ν₀ = experimental_data.ν[11] # based on Kleinheerenbrink defintion of when perching starts
    χ₀ = experimental_data.χ[11]
    ζ₀ = experimental_data.ζ[11]
    χ₁ = experimental_data.χ[end]
    ζ₁ = experimental_data.ζ[end]
    return @def begin
        T ∈ R, variable
        t ∈ [0, T], time
        # ν̇ = -CD(α(t), p) * ν(t)^2 - sin(γ(t)) + CT(t)
        # γ̇ = CL(α(t), p) * ν(t) - cos(γ(t)) / ν(t)
        x = (ν, γ, χ, ζ) ∈ R⁴, state
        u = (α, CT) ∈ R², control
        ∂(ν)(t) == -CD(α(t), p) * ν(t)^2 - sin(γ(t)) + CT(t)
        ∂(γ)(t) == CL(α(t), p) * ν(t) - cos(γ(t)) / ν(t)
        ∂(χ)(t) == ν(t) * cos(γ(t))
        ∂(ζ)(t) == ν(t) * sin(γ(t))
        # ẋ(t) == [ν̇, γ̇, ν(t) * cos(γ(t)), ν(t) * sin(γ(t))]
        ν(t) ≥ 0.01
        -π / 2 ≤ γ(t) ≤ π / 2
        0.0 ≤ α(t) ≤ π / 2
        CTmin ≤ CT(t) ≤ CTmax
        χ(0) == χ₀
        ζ(0) == ζ₀
        ν(0) == ν₀
        lb(χ₁, pos_bounds) ≤ χ(T) ≤ ub(χ₁, pos_bounds)
        lb(ζ₁, pos_bounds) ≤ ζ(T) ≤ ub(ζ₁, pos_bounds)
        (λ * ν(T)) + (1 - λ) * ∫(post_stall_x(ν(t), α(t), stall_angle)) → min
    end
end

"""
dc_weighted_min_speed_poststall_time(p, λ; CTmin=0.0, CTmax=1.0, pos_bounds=0.05, stall_angle=deg2rad(25))

Create optimal control problem which minimises weighted cost between touchdown speed and post-stall *time*. Weighted by parameter `λ`

`pos_bounds` chooses the tolerance to the measured touchdown position from the experimental data.
"""

function dc_weighted_min_speed_poststall_time(p, λ; CTmin=0.0, CTmax=1.0, pos_bounds=0.05, stall_angle=deg2rad(25), νmin=0.0, νmax=1.0, Tmin=0.0, Tmax=1.0)
    ν₀ = experimental_data.ν[11] # based on Kleinheerenbrink defintion of when perching starts
    χ₀ = experimental_data.χ[11]
    ζ₀ = experimental_data.ζ[11]
    χ₁ = experimental_data.χ[end]
    ζ₁ = experimental_data.ζ[end]
    return @def begin
        T ∈ R, variable
        t ∈ [0, T], time
        x = (ν, γ, χ, ζ) ∈ R⁴, state
        u = (α, CT) ∈ R², control
        ∂(ν)(t) == -CD(α(t), p) * ν(t)^2 - sin(γ(t)) + CT(t)
        ∂(γ)(t) == CL(α(t), p) * ν(t) - cos(γ(t)) / ν(t)
        ∂(χ)(t) == ν(t) * cos(γ(t))
        ∂(ζ)(t) == ν(t) * sin(γ(t))
        # ẋ(t) == [ν̇, γ̇, ν(t) * cos(γ(t)), ν(t) * sin(γ(t))]
        ν(t) ≥ 0.01
        -π / 2 ≤ γ(t) ≤ π / 2
        0.0 ≤ α(t) ≤ π / 2
        CTmin ≤ CT(t) ≤ CTmax
        χ(0) == χ₀
        ζ(0) == ζ₀
        ν(0) == ν₀
        lb(χ₁, pos_bounds) ≤ χ(T) ≤ ub(χ₁, pos_bounds)
        lb(ζ₁, pos_bounds) ≤ ζ(T) ≤ ub(ζ₁, pos_bounds)
        term_cost = (ν(T) - νmin) / (νmax - νmin)
        intg_cost_constant = ((1 - λ) * Tmin) / (Tmax - Tmin)
        mayer = (λ * term_cost - intg_cost_constant)
        mayer + ((1 - λ) / (Tmax - Tmin)) * ∫(post_stall_x(1.0, α(t), stall_angle)) → min
    end
end
