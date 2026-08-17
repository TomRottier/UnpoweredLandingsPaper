# convert solution to dimensional form
nominal_speed(m, S, ρ, g) = √(2m * g / (ρ * S))
nominal_time(m, S, ρ, g) = nominal_speed(m, S, ρ, g) / g
nominal_length(m, S, ρ) = 2m / (ρ * S)
nominal_speed(mdl::DimensionalModel) = nominal_speed(mdl.m, mdl.S, mdl.ρ, mdl.g)
nominal_time(mdl::DimensionalModel) = nominal_time(mdl.m, mdl.S, mdl.ρ, mdl.g)
nominal_length(mdl::DimensionalModel) = nominal_length(mdl.m, mdl.S, mdl.ρ)

struct DimensionalSolution
    t::Vector{Float64}
    v::Vector{Float64}
    γ::Vector{Float64}
    x::Vector{Float64}
    z::Vector{Float64}
end

function dimensionalise_sol(model::DimensionalModel, sol)
    # nominal values
    t̂ = nominal_time(model)
    v̂ = nominal_speed(model)
    l̂ = nominal_length(model)

    # dimensionalise solution
    t = t̂ * sol.t
    v = v̂ * sol[1, :]
    γ = sol[2, :]
    x = l̂ * sol[3, :]
    z = l̂ * sol[4, :]

    return Solution(t, v, γ, x, z)
end

