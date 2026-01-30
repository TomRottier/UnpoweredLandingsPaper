module GlideLandingPaper

using CairoMakie, PointMass, OrdinaryDiffEqTsit5, NonlinearSolve, DelimitedFiles, Folds
using QuadGK: quadgk
using OptimalControl: OptimalControl, @def, solve, time_grid, state, control, variable, objective, CTModels.success
using NLPModelsIpopt, MadNLP, MadNLPMumps
using LinearAlgebra: ⋅, dot
using CairoMakie.GeometryBasics: Polygon
import CairoMakie.Makie.SpecApi as S
using ColorSchemes: colorschemes
using Logging

include("constants.jl")
include("dynamics.jl")
include("control.jl")
include("dc.jl")
include("plotting.jl")

include("../figs/methods.jl")
include("../figs/dynamics.jl")
include("../figs/optcontrol.jl")
include("../figs/constraints.jl")
include("../figs/perching.jl")
include("../figs/sensitivity.jl")


# set plot theme globally at init
__init__() = set_theme!(plottheme)

const FIG_SAVE_PATH = joinpath(dirname(dirname(@__DIR__)), "figs")

"""
    create_figures(figs; kwargs...)

Generates then saves figures specified in `figs`.

`figs` is a dictionary with keys being the figure name and value a function to generate the figure, `ft`, specifies the filetype to save as, `kwargs` are passed on to `Makie.save`
"""
function create_figures(figs; ft = "pdf", kwargs...)
    for (name, fcn) in figs
        @info "generating $name"
        with_logger(ConsoleLogger(Error)) do
            save(joinpath(FIG_SAVE_PATH, string(name, ".", ft)), fcn(); kwargs...)
        end
    end
    return
end


plotfcns = [
    freebodydiagram,
    dynamics_phasespace,
    optcontrol,
    constraints,
    perching,
    sensitivity,
]
figs_dict = Dict(string(fcn) => fcn for fcn in plotfcns)

function (@main)(; ft = "pdf", kwargs...)
    return create_figures(figs_dict; ft, kwargs)
end


export create_figures, figs_dict

#### precompile
# using PrecompileTools: @compile_workload

# @compile_workload begin
#     foreach(Base.Fix2(invoke, Tuple{}), plotfcns)
# end

end # module GlideLandingPaper
