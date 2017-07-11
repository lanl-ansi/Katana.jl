module Katana

using MathProgBase
using ConicNonlinearBridge
using JuMP

include("separators.jl")
include("algorithms.jl")

# parameter passing between KatanaSolver and KatanaNonlinearModel
immutable KatanaModelParams
    f_tol        :: Float64 # feasibility tolerance
    aux_lb       :: Float64 # auxiliary variable lower bound
    aux_ub       :: Float64 # '' upper bound
    iter_cap     :: Int64   # iteration cap
    separator    :: AbstractKatanaSeparator # separation oracle
end

include("solver.jl")
include("model.jl")
include("util.jl")

end
