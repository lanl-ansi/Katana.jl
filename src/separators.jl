export AbstractKatanaSeparator, KatanaFirstOrderSeparator, initialize!, precompute!, gencut!

abstract type AbstractKatanaSeparator end

"""
    initialize!(sep::AbstractKatanaSeparator, linear_model, num_constr, oracle::MathProgBase.AbstractNLPEvaluator)

Initialise an instance of a subtype of `AbstractKatanaSeparator`. This method is called by `loadproblem!` and MUST be overridden.
`linear_model` is the internal linear model of the `KatanaNonlinearModel`. `num_constr` is the number of constraints in the problem,
excluding the epigraph constraint on the objective function (see solver implementation description). `oracle` can be queried
for first- and second- derivative information and must be initialised in this method (see MathProgBase documentation on nonlinear models).
"""
initialize!(sep::AbstractKatanaSeparator, linear_model, num_constr, oracle::MathProgBase.AbstractNLPEvaluator) = error("Not implemented: Katana.initialize!")

"""
    gencut!(separator::AbstractKatanaSeparator, xstar, i)

Generate a cut given an LP solution `xstar` for constraint `i`. If `i` is greater than `num_constr` passed in `initialize!`,
treat the objective as a constraint and generate a cut for it. This method MUST be overridden for a subtype of `AbstractKatanaSeparator`.
"""
gencut!(separator::AbstractKatanaSeparator, xstar, i) = error("Not implemented: Katana.gencut!")

"""
    precompute!(separator :: AbstractKatanaSeparator, xstar)

Implement this method for a subtype of AbstractKatanaSeparator if your separator might only need to evaluate
certain information once for all constraints using a solution from the internal LP model.

`xstar` is the solution vector from the `KatanaNonlinearModel`'s LP model, including the added auxiliary variable
representing the value of the objective function (which is itself evaluated as a constraint).
"""
precompute!(separator :: AbstractKatanaSeparator) = nothing

"""
An implementation of `AbstractKatanaSeparator` for any first-order cutting algorithm.
"""
type KatanaFirstOrderSeparator <: AbstractKatanaSeparator
    sp_cols :: Vector{Vector{Int}} # column indices in Jacobian by row
    sp_col_inds :: Vector{Vector{Int}} # original indices in sparse Jacobian vector

    linear_model :: JuMP.Model # JuMP LP model used by KatanaNonlinearModel
    oracle :: MathProgBase.AbstractNLPEvaluator # evaluator queried as oracle

    jac :: Vector{Float64} # jacobian sparse matrix
    g   :: Vector{Float64} # vector of constraint values returned by MathProgBase.eval_g()

    algo # cutting-plane generator that takes: separator, point and constraint index and returns an AffExpr

    function KatanaFirstOrderSeparator(algo) = (s = new(); s.algo = algo) # TODO default algo?
end

# Implements initialize! for KatanaFirstOrderSeparators. This initialises the oracle with necessary features
#  and processes the sparsity structure of the Jacobian to allow faster lookups of non-zero entries by row.
function initialize!(sep          :: KatanaFirstOrderSeparator,
                     linear_model :: JuMP.Model,
                     num_constr   :: Int,
                     oracle       :: MathProgBase.AbstractNLPEvaluator)
    sep.linear_model = linear_model

    MathProgBase.initialize(oracle, [:Grad, :Jac, :JacVec]) # certain cut algos may need jacobian-vector products?
    sep.oracle = oracle

    # map info of Jacobian J_g to allow faster lookups of nonzero entries by row:
    sp_rows, sp_cols = MathProgBase.jac_structure(oracle)
    sep.sp_cols = Vector{Int}[ [] for i=1:num_constr] # map a row to vector of nonzero columns' indices
    N = length(sp_rows) # number of sparse entries
    for ind in 1:N
        i,j = sp_rows[ind],sp_cols[ind]
        push!(sep.sp_cols[i], j) # column j is nonzero in row i
        push!(sep.sp_col_inds[i], ind)
    end

    sep.jac = zeros(N)
    sep.g = zeros(num_constr)
end

# Implements precompute! for KatanaFirstOrderSeparators. This evaluates every constraint at the point xstar and computes
#  and stores the sparse Jacobian at xstar
function precompute!(separator :: KatanaFirstOrderSeparator, xstar)
    MathProgBase.eval_jac_g(separator.oracle , separator.jac, xstar[1:end-1]) # hopefully variable ordering is consistent with MPB
    MathProgBase.eval_g(separator.oracle, separator.g, xstar[1:end-1]) # evaluate constraints
end

function gencut!(separator :: KatanaFirstOrderSeparator, xstar, i)
    cut = separator.algo(separator, xstar, i) # TODO
    # TODO add cut to LP model
end
