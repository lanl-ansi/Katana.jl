export KatanaSolver

"""
A solver for convex NLPs that uses cutting-planes to approximate a convex feasible set.
"""
type KatanaSolver <: MathProgBase.AbstractMathProgSolver
    lp_solver    :: MathProgBase.AbstractMathProgSolver
    features     :: Vector{Symbol}
    model_params :: KatanaModelParams
end

"""
    KatanaSolver(lp_solver :: MathProgBase.AbstractMathProgSolver;
                 separator = KatanaFirstOrderSeparator(),
                 features = Vector{Symbol}(),
                 f_tol :: Float64 = 1e-6,
                 cut_coef_rng :: Float64 = 1e9,
                 log_level :: Int = 10,
                 iter_cap :: Int = 10000,
                 obj_eps :: Float64 = -1.0)

Construct a `KatanaSolver` with feasibility tolerance `f_tol`, a maximum coefficient range per cut `cut_coef_rng` and an
iteration cap specifying the maximum number of rounds of LP solves + cut generation in `iter_cap`. Print out solver progress
every `log_level` number of iterations, or suppress output with `log_level=0`. Use `obj_eps` as a stopping criterion: when the objective function evaluated by the LP has changed by less than `obj_eps` relative to the previous iteration.

`cut_coef_rng` is used to round-off close-to-zero coefficients in generated cuts.

The `separator` is any implementing subtype of `AbstractKatanaSeparator`. It serves as the separation oracle
used by the solver. The default is a first order separator that generates a single Newton cut per constraint.

The `features` vector is a list of optional features to enable in the solver. Currently supported are
* `:VisData` ``-`` Internal model logs actions to be exported and visualised
"""
function KatanaSolver(lp_solver::MathProgBase.AbstractMathProgSolver;
                      separator = KatanaFirstOrderSeparator(),
                      features = Vector{Symbol}(),
                      f_tol    :: Float64 = 1e-6,
                      cut_coef_rng :: Float64 = 1e9,
                      log_level :: Int = 10,
                      iter_cap :: Int     = 10000,
                      obj_eps  :: Float64 = -1.0)
    return KatanaSolver(lp_solver, features, KatanaModelParams(f_tol, iter_cap, log_level, cut_coef_rng, obj_eps, separator))
end

# this bridge should make lp/qp models act like nlp models
MathProgBase.LinearQuadraticModel(s::KatanaSolver) = MathProgBase.NonlinearToLPQPBridge(MathProgBase.NonlinearModel(s))
