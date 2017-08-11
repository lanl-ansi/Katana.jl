module Katana

using MathProgBase
using ConicNonlinearBridge
using JuMP

include("nlpeval.jl")
include("separators.jl")
include("algorithms.jl")

# parameter passing between KatanaSolver and KatanaNonlinearModel
immutable KatanaModelParams
    f_tol        :: Float64 # feasibility tolerance
    iter_cap     :: Int   # iteration cap
    presolve_cap :: Int     # cap on number of times to run bounding routine in presolve step
    log_level    :: Int     # printout frequency (every iter_freq iterations)
    cut_coef_rng :: Float64 # max coefficient range per cut
    obj_eps      :: Float64 # stopping criteria on objective delta
    separator    :: AbstractKatanaSeparator # separation oracle
end

include("solver.jl")
include("model.jl")
include("util.jl")

end
