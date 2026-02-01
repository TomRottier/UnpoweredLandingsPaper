module PointMass

export Parameters, Aerodynamics, Control, DimensionalModel, Solution # types
export ν̇, γ̇, χ̇, ζ̇, dynamics!, dynamics, full_dynamics!, full_dynamics # dynamics
export kinetic_energy, potential_energy, total_energy, aerodynamic_power, gravitational_power, thrust_power
export TrigAerodynamics, LinearAerodynamics, LinearSmoothAerodynamics # aerodynamic models
export CL, CD, inverse_CL, inverse_CD, alpha_CL_max, alpha_CD_max, alpha_CL_CD_max, ∂CL_∂α, ∂CD_∂α, min_CD, max_CD, max_CL # aerodynamic functions
export alpha, thrust, thrust_angle, ConstantAlpha, LinearAlpha, InterpolatedControl # control
export nominal_time, nominal_speed, nominal_length, dimensionalise_sol # dimensional conversion functions

include("types.jl")
include("aerodynamics.jl")
include("control.jl")
include("dynamics.jl")
include("functions.jl")
include("analysis.jl")

end # module PointMass
