CL(α, p::Parameters{<:Aerodynamics,<:Control}) = CL(α, p.aerodynamics)
CD(α, p::Parameters{<:Aerodynamics,<:Control}) = CD(α, p.aerodynamics)
alpha_CL_max(p::Parameters{<:Aerodynamics,<:Control}) = alpha_CL_max(p.aerodynamics)
alpha_CD_max(p::Parameters{<:Aerodynamics,<:Control}) = alpha_CD_max(p.aerodynamics)
alpha_CL_CD_max(p::Parameters{<:Aerodynamics,<:Control}) = alpha_CL_CD_max(p.aerodynamics)
# inverse_CL(cl, p::Parameters{<:Aerodynamics,<:Control}) = inverse_CL(cl, p.aerodynamics)
inverse_CD(cd, p::Parameters{<:Aerodynamics,<:Control}) = inverse_CD(cd, p.aerodynamics)
∂CL_∂α(α, p::Parameters{<:Aerodynamics,<:Control}) = ∂CL_∂α(α, p.aerodynamics)
∂CD_∂α(α, p::Parameters{<:Aerodynamics,<:Control}) = ∂CD_∂α(α, p.aerodynamics)
min_CD(p::Parameters{<:Aerodynamics,<:Control}) = min_CD(p.aerodynamics)
max_CD(p::Parameters{<:Aerodynamics,<:Control}) = max_CD(p.aerodynamics)
max_CL(p::Parameters{<:Aerodynamics,<:Control}) = max_CL(p.aerodynamics)

## fallback case
CL(α, p) = p.CLmax * sin(2α)
CD(α, p) = p.CD₀ + p.CDmax * sin(α)^2

## simple trigonometric functions for CL,CD
struct TrigAerodynamics <: Aerodynamics
    CLmax::Float64
    CD₀::Float64
    CDmax::Float64
end
TrigAerodynamics() = TrigAerodynamics(1.5, 0.03, 2.0)
CL(α, p::TrigAerodynamics) = p.CLmax * sin(2α)
CD(α, p::TrigAerodynamics) = p.CD₀ + p.CDmax * sin(α)^2
alpha_CL_max(p::TrigAerodynamics) = π / 4
alpha_CD_max(p::TrigAerodynamics) = π / 2
alpha_CL_CD_max(p::TrigAerodynamics) = atan(sqrt(p.CD₀ / (p.CD₀ + p.CDmax)))
inverse_CD(cd, p::TrigAerodynamics) = asin(sqrt((cd - p.CD₀) / p.CDmax))
# inverse_CL(cl, p::TrigAerodynamics) = cl ≤ p.CLmax ? 0.5asin(cl / p.CLmax) : π / 2 - 0.5asin(cl / p.CLmax)
∂CL_∂α(α, p::TrigAerodynamics) = 2p.CLmax * cos(2α)
∂CD_∂α(α, p::TrigAerodynamics) = p.CDmax * sin(2α)
min_CD(p::TrigAerodynamics) = p.CD₀
max_CD(p::TrigAerodynamics) = p.CD₀ + p.CDmax
max_CL(p::TrigAerodynamics) = p.CLmax

## include linear pre-stall effects for lift and quadratic for drag
struct LinearAerodynamics <: Aerodynamics
    slope::Float64
    CD₀::Float64
    k::Float64
    αstall::Float64
    poststall::TrigAerodynamics

    LinearAerodynamics(slope, CD₀, k, αstall, CLmax, CDmax) = new(slope, CD₀, k, αstall, TrigAerodynamics(CLmax, CD₀, CDmax))
end
LinearAerodynamics(slope, k, αstall, post::TrigAerodynamics) = LinearAerodynamics(slope, post.CD₀, k, αstall, post.CLmax, post.CDmax)
LinearAerodynamics() = LinearAerodynamics(2π, 0.03, deg2rad(15), TrigAerodynamics(1.5, 0.03, 2.0))
CL(α, p::LinearAerodynamics) = α ≤ p.αstall ? p.slope * α : CL(α, p.poststall)
CD(α, p::LinearAerodynamics) = α ≤ p.αstall ? p.CD₀ + p.k * p.slope^2 * α^2 : CD(α, p.poststall)
alpha_CL_max(p::LinearAerodynamics) = p.slope * p.αstall > p.CLmax ? p.αstall : π / 4
alpha_CD_max(p::LinearAerodynamics) = p.k * p.slope^2 * p.αstall^2 > p.CDmax ? p.αstall : π / 2
alpha_CL_CD_max(p::LinearAerodynamics) = p.k * p.slope^2 > p.CDmax ? sqrt(p.CD₀ / p.k * p.a^2) : atan(sqrt(p.CD₀ / (p.CD₀ + p.CDmax)))
∂CL_∂α(α, p::LinearAerodynamics) = α < p.αstall ? p.slope : ∂CL_∂α(α, p.poststall)
∂CD_∂α(α, p::LinearAerodynamics) = α < p.αstall ? 2p.k * p.slope^2 * α : ∂CD_∂α(α, p.poststall)
min_CD(p::LinearAerodynamics) = p.CD₀
max_CD(p::LinearAerodynamics) = p.CD₀ + p.poststall.CDmax
max_CL(p::LinearAerodynamics) = p.poststall.CLmax
inverse_CD(cd, p::LinearAerodynamics) = cd ≤ CD(p.αstall, p) ? sqrt((cd - p.CD₀) / (p.k * p.slope^2)) : inverse_CD(cd, p.poststall)
# inverse_CL(cl, p::LinearAerodynamics) = cl ≤ CL(p.αstall, p) ? cl / p.slope : inverse_CL(cl, p.poststall)

## smoothed version of linear pre-stall, taken from Li et al., 2022 J. Fluid. Mech.
struct LinearSmoothAerodynamics <: Aerodynamics
    CL₁::Float64 # approx slope of pre-stall CL
    CL₂::Float64 # max post-stall CL
    CD₀::Float64 # skin friction drag
    CD₁::Float64 # approx slope of pre-stall CD
    CD₂::Float64 # max post stall CD
    αstall::Float64 # α₀ in Li et al.
    δ::Float64 # smoothness of transition between pre-stall and stall
end
LinearSmoothAerodynamics() = LinearSmoothAerodynamics(5.3, 0.95, 0.1, 5.0, 1.9, deg2rad(14), deg2rad(6)) # values from Li et al
LinearSmoothAerodynamics(p::TrigAerodynamics) = LinearSmoothAerodynamics(2π, p.CLmax, p.CD₀, 0.1, p.CDmax, deg2rad(14), deg2rad(6))
f̃(α, αstall, δ) = 0.5 * (1 - tanh((α - αstall) / δ)) # α ∈ [0, π/2]
C̃L(α, CL₁, CL₂, αstall, δ) = f̃(α, αstall, δ) * CL₁ * sin(α) + (1 - f̃(α, αstall, δ)) * CL₂ * sin(2α) # α ∈ [0, π/2]
C̃D(α, CD₀, CD₁, CD₂, αstall, δ) = f̃(α, αstall, δ) * (CD₀ + CD₁ * sin(α)^2) + (1 - f̃(α, αstall, δ)) * (CD₂ * sin(α)^2) # α ∈ [0, π/2]
∂f̃_∂α(α, αstall, δ) = -sech((αstall - α) / δ)^2 / 2δ
# CL(α, p::LinearSmoothAerodynamics) = begin
#     if -π ≤ α ≤ -π / 2
#         C̃L(π - abs(α), p.CL₁, p.CL₂, p.αstall, p.δ)
#     elseif -π / 2 ≤ α ≤ 0.0
#         -C̃L(abs(α), p.CL₁, p.CL₂, p.αstall, p.δ)
#     elseif 0.0 ≤ α ≤ π / 2
#         C̃L(α, p.CL₁, p.CL₂, p.αstall, p.δ)
#     elseif π / 2 ≤ α ≤ π
#         -C̃L(π - α, p.CL₁, p.CL₂, p.αstall, p.δ)
#     else
#         error("α out of bounds")
#     end
# end
# CD(α, p::LinearSmoothAerodynamics) = begin
#     if -π ≤ α ≤ -π / 2
#         C̃D(π - abs(α), p.CD₀, p.CD₁, p.CD₂, p.αstall, p.δ)
#     elseif -π / 2 ≤ α ≤ 0.0
#         C̃D(abs(α), p.CD₀, p.CD₁, p.CD₂, p.αstall, p.δ)
#     elseif 0.0 ≤ α ≤ π / 2
#         C̃D(α, p.CD₀, p.CD₁, p.CD₂, p.αstall, p.δ)
#     elseif π / 2 ≤ α ≤ π
#         C̃D(π - α, p.CD₀, p.CD₁, p.CD₂, p.αstall, p.δ)
#     else
#         error("α out of bounds")
#     end
# end
CL(α, p::LinearSmoothAerodynamics) = C̃L(α, p.CL₁, p.CL₂, p.αstall, p.δ)
CD(α, p::LinearSmoothAerodynamics) = C̃D(α, p.CD₀, p.CD₁, p.CD₂, p.αstall, p.δ)
min_CD(p::LinearSmoothAerodynamics) = p.CD₀
max_CD(p::LinearSmoothAerodynamics) = p.CD₂
max_CL(p::LinearSmoothAerodynamics) = max(CL(p.αstall, p), CL(π / 4, p))
∂CL_∂α(α, p::LinearSmoothAerodynamics) =
    ∂f̃_∂α(α, p.αstall, p.δ) * (p.CL₁ * sin(α) - p.CL₂ * sin(2α)) + f̃(α, p.αstall, p.δ) * (p.CL₁ * cos(α) - 2p.CL₂ * cos(2α)) + 2p.CL₂ * cos(2α)
∂CD_∂α(α, p::LinearSmoothAerodynamics) =
    ∂f̃_∂α(α, p.αstall, p.δ) * (p.CD₀ + (p.CD₁ - p.CD₂) * sin(α)^2) + f̃(α, p.αstall, p.δ) * sin(2α) * (p.CD₁ - p.CD₂) + p.CD₂ * sin(2α)

