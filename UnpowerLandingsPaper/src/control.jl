export OptControl

"""
    OptControl <: PointMass.Control

Optimal control strategy that reduces speed as much as possible by minimising ratio of γ̇/ν̇.

"""
struct OptControl <: Control end

# the numerator of ∂γ̇/∂ν̇
∂γ̇_∂ν̇(α, ν, γ, p) =
    ∂CL_∂α(α, p) * (-CD(α, p) * ν^2 - sin(γ)) +
    ∂CD_∂α(α, p) * ν * (CL(α, p) * ν - cos(γ) / ν)

function isoutofdomain(u, p::Parameters{<:Aerodynamics,OptControl}, t)
    ν, γ = u
    cd0 = -sin(γ) / ν^2
    cd0 ≈ max_CD(p) && return false

    if cd0 < min_CD(p)
        return false
    elseif min_CD(p) ≤ cd0 ≤ max_CD(p)
        αmin = inverse_CD(min(cd0, max_CD(p)), p)
        γ̇ = CL(αmin, p) * ν - cos(γ) / ν
        if γ̇ > 0.0
            return true
        else
            return false
        end
    else
        return true
    end
end

function PointMass.alpha(u, p::Parameters{<:Aerodynamics,OptControl}, t)
    if !isoutofdomain(u, p, t)
        ν, γ = u
        cd0 = -sin(γ) / ν^2
        if cd0 < min_CD(p) # ν̇ < 0
            αmin = zero(eltype(u))
        elseif min_CD(p) ≤ cd0 ≤ max_CD(p) || cd0 ≈ max_CD(p) # ν̇ can = 0
            αmin = inverse_CD(min(cd0, max_CD(p)), p) # only consider those α which cause ν̇ ≤ 0
            αmin ≈ π / 2 && return eltype(u)(π / 2)
        end
        prob = IntervalNonlinearProblem((u, p) -> ∂γ̇_∂ν̇(u, ν, γ, p), (αmin, π / 2), p)
        sol = solve(prob, ITP())
        return sol.u
    else
        return eltype(u)(10)
    end
end
PointMass.thrust(u, p::OptControl, t) = 0.0
PointMass.thrust_angle(u, p::OptControl, t) = 0.0

function acceleration_nullcline(ν, p::Parameters{<:Aerodynamics,OptControl})
    ff(u, p) = begin
        γ = u
        α = alpha([ν, γ], p, nothing)
        α == 10 && return Inf
        return -CD(α, p) * ν^2 - sin(γ)
    end
    1.0 + ff(-π / 2, p) ≈ 1.0 && return -π / 2 # using ≈ 0.0 is equivalent to == 0 (see isapprox docs)
    prob = IntervalNonlinearProblem(ff, (-π / 2, 0.1), p)
    sol = solve(prob, ITP(), abstol=1.0e-6, reltol=1.0e-6)
    return sol.u
end

function acceleration_nullcline(p::Parameters{<:Aerodynamics,OptControl}; n)
    νs = [
        range(0.0, sqrt(1 / max_CD(p)), n ÷ 2);
        range(sqrt(1 / max_CD(p)), sqrt(1 / min_CD(p)), n ÷ 2)
    ] # gives better coverage
    γs = map(Base.Fix2(acceleration_nullcline, p), νs)
    return (νs, γs)
end

function flightpath_nullcline(ν, p)
    ff(u, p) = begin
        γ = u
        α = alpha([ν, γ], p, nothing)
        α == 10 && return Inf
        return CL(α, p) * ν - cos(γ) / ν
    end
    prob = IntervalNonlinearProblem(ff, (0.0, 1.6), p)
    sol = solve(prob, ITP(), abstol=1.0e-12, reltol=1.0e-12)
    return sol.u
end

function flightpath_nullcline(p::Parameters{<:Aerodynamics,OptControl}; n)
    νs = range(1.0e-8, sqrt(1 / max_CL(p)), length=n ÷ 2)
    γs = map(Base.Fix2(flightpath_nullcline, p), νs)
    return ([νs; reverse(νs)], [γs; reverse(γs)])
end


"""
    PreStallAerodynamics <: Aerodynamics


Same as LinearSmootherAerodynamics from PointMass

Fields:
    `CL₁`: pre-stall CL slope
    `CL₂`: max post-stall CL
    `CD₀`: skin friction drag
    `CD₁`: approx slope of pre-stall CD
    `CD₂`: max post stall CD
    `αstall`: α₀ in Li et al.
    `δ`: smoothness of transition between pre-stall and stall

"""
struct PreStallAerodynamics <: Aerodynamics
    CL₁::Float64 # approx slope of pre-stall CL
    CL₂::Float64 # max post-stall CL
    CD₀::Float64 # skin friction drag
    CD₁::Float64 # approx slope of pre-stall CD
    CD₂::Float64 # max post stall CD
    αstall::Float64 # α₀ in Li et al.
    δ::Float64 # smoothness of transition between pre-stall and stall
end

PreStallAerodynamics() =
    PreStallAerodynamics(2π, 1.5, 0.03, 0.1, 2.0, deg2rad(20), deg2rad(6))
PreStallAerodynamics(p::TrigAerodynamics) =
    (2π, p.CLmax, p.CD₀, 0.1, p.CDmax, deg2rad(20), deg2rad(6))


f̃(α, αstall, δ) = 0.5 * (1 - tanh((α - αstall) / δ)) # α ∈ [0, π/2]
C̃L(α, CL₁, CL₂, αstall, δ) =
    f̃(α, αstall, δ) * CL₁ * sin(α) + (1 - f̃(α, αstall, δ)) * CL₂ * sin(2α) # α ∈ [0, π/2]
C̃D(α, CD₀, CD₁, CD₂, αstall, δ) =
    f̃(α, αstall, δ) * (CD₀ + CD₁ * sin(α)^2) + (1 - f̃(α, αstall, δ)) * (CD₂ * sin(α)^2) # α ∈ [0, π/2]
∂f̃_∂α(α, αstall, δ) = -sech((αstall - α) / δ)^2 / 2δ
PointMass.CL(α, p::PreStallAerodynamics) = C̃L(α, p.CL₁, p.CL₂, p.αstall, p.δ)
PointMass.CD(α, p::PreStallAerodynamics) = C̃D(α, p.CD₀, p.CD₁, p.CD₂, p.αstall, p.δ)
PointMass.min_CD(p::PreStallAerodynamics) = p.CD₀
PointMass.max_CD(p::PreStallAerodynamics) = p.CD₂
PointMass.max_CL(p::PreStallAerodynamics) = max(CL(p.αstall, p), CL(π / 4, p))
PointMass.∂CL_∂α(α, p::PreStallAerodynamics) =
    ∂f̃_∂α(α, p.αstall, p.δ) * (p.CL₁ * sin(α) - p.CL₂ * sin(2α)) +
    f̃(α, p.αstall, p.δ) * (p.CL₁ * cos(α) - 2p.CL₂ * cos(2α)) +
    2p.CL₂ * cos(2α)
PointMass.∂CD_∂α(α, p::PreStallAerodynamics) =
    ∂f̃_∂α(α, p.αstall, p.δ) * (p.CD₀ + (p.CD₁ - p.CD₂) * sin(α)^2) +
    f̃(α, p.αstall, p.δ) * sin(2α) * (p.CD₁ - p.CD₂) +
    p.CD₂ * sin(2α)

PointMass.inverse_CD(cd, p::PreStallAerodynamics) =
    solve(IntervalNonlinearProblem((α, p) -> CD(α, p) - cd, (0.0, π / 2), p), ITP()).u
