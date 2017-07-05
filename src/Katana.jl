module Katana

using MathProgBase
using ConicNonlinearBridge
using JuMP

export KatanaSolver

# parameter passing between KatanaSolver and KatanaNonlinearModel
immutable KatanaModelParams
    f_tol        :: Float64 # feasibility tolerance
    aux_lb       :: Float64 # auxiliary variable lower bound
    aux_ub       :: Float64 # '' upper bound
    iter_cap     :: Int64   # iteration cap
end

include("solver.jl")
include("model.jl")

end
