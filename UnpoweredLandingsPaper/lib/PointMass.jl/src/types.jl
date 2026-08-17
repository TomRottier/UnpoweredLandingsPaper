# abstract type for aerodynamics models
abstract type Aerodynamics end

# abstract type for all control models
abstract type Control end

# struct for parameters for simulation
struct Parameters{T1<:Aerodynamics,T2<:Control}
    aerodynamics::T1
    control::T2
end
Parameters() = Parameters(TrigAerodynamics(1.5, 0.03, 2.0), ConstantControlNoThrust(0.2))

# struct to hold parameters for dimensional model
struct DimensionalModel{T}
    m::T
    S::T
    ρ::T
    g::T
end

# default parameters
DimensionalModel(m, S, ρ=1.225, g=9.81) = DimensionalModel(m, S, ρ, g)
DimensionalModel() = DimensionalModel(0.3, 0.04, 1.225, 9.81)

# own solution type
struct Solution{T}
    t::Vector{T}
    v::Vector{T}
    γ::Vector{T}
    x::Vector{T}
    z::Vector{T}
end