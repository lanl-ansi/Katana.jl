export KatanaSolver

"""
docs go here
"""
type KatanaSolver <: MathProgBase.AbstractMathProgSolver
    lp_solver    :: MathProgBase.AbstractMathProgSolver
    features     :: Vector{Symbol}
    model_params :: KatanaModelParams
end

"""
    KatanaSolver(lp_solver::MathProgBase.AbstractMathProgSolver;
                 features :: Vector{Symbol} = [],
                 f_tol    :: Float64 = 1e-6,
                 aux_lb   :: Float64 = -1e6,
                 aux_ub   :: Float64 = 1e6,
                 iter_cap :: Int     = 10000)

Construct a `KatanaSolver` with feasibility tolerance `f_tol`, lower and upper
bounds on auxiliary variables `aux_lb` and `aux_ub`, and an iteration cap specifying the number
of rounds of LP solves + cut generation in `iter_cap`.

The `features` vector is a list of optional features to enable in the solver. Currently supported are
* `:VisData` ``-`` Internal model logs actions to be exported and visualised
"""
function KatanaSolver(lp_solver::MathProgBase.AbstractMathProgSolver;
                      features = Vector{Symbol}(),
                      f_tol    :: Float64 = 1e-6,
                      aux_lb   :: Float64 = -1e6,
                      aux_ub   :: Float64 = 1e6,
                      iter_cap :: Int     = 10000)
    return KatanaSolver(lp_solver, features, KatanaModelParams(f_tol, aux_lb, aux_ub, iter_cap))
end

# this bridge should make lp/qp models act like nlp models
MathProgBase.LinearQuadraticModel(s::KatanaSolver) = MathProgBase.NonlinearToLPQPBridge(MathProgBase.NonlinearModel(s))
