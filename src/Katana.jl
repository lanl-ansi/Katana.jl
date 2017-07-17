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
    iter_cap     :: Int64   # iteration cap
    log_level    :: Int     # printout frequency (every iter_freq iterations)
    cut_coef_rng :: Float64 # max coefficient range per cut
    separator    :: AbstractKatanaSeparator # separation oracle
end

include("solver.jl")
include("model.jl")
include("util.jl")

end
