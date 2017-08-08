export AbstractKatanaSeparator, KatanaFirstOrderSeparator, KatanaProjectionSeparator, initialize!, precompute!, gencut!

using Compat

@compat abstract type AbstractKatanaSeparator end

"""
    initialize!(sep::AbstractKatanaSeparator, linear_model, num_var, num_constr, oracle::MathProgBase.AbstractNLPEvaluator)

Initialise an instance of a subtype of `AbstractKatanaSeparator`. This method is called by `loadproblem!` and MUST be overridden.
`linear_model` is the internal linear model of the `KatanaNonlinearModel`.

`num_var` is the number of solution variables, as passed by Katana. `num_constr` is the number of
constraints in the problem, as passed by Katana. See solver implementation documentation for details.

`oracle` can be queried for first- and second- derivative information and must be initialised in this method (see MathProgBase
documentation on nonlinear models).
"""
initialize!(sep::AbstractKatanaSeparator, linear_model, num_var, num_constr, oracle::MathProgBase.AbstractNLPEvaluator) = error("Not implemented: Katana.initialize!")

"""
    gencut!(sep::AbstractKatanaSeparator, xstar, i)

Generate a list of cuts given some input `a` for constraint `i`.
This method MUST be overridden for a subtype of `AbstractKatanaSeparator`.

Return a JuMP.AffExpr object.
"""
gencuts(sep::AbstractKatanaSeparator, a, i) = error("Not implemented: Katana.gencuts!")

# compatibility
gencut(sep::AbstractKatanaSeparator, a, i) = (cs = gencuts(sep,a,i); cs[1])

"""
    isconstrsat(sep::AbstractKatanaSeparator, i, lb, ub, f_tol)

Returns true if the ith constraint is satisfied for the given bounds and tolerance. This method MUST be overriden for a subtype of
`AbstractKatanaSeparator` as querying evaluated constraints is implementation-dependent.
"""
isconstrsat(sep::AbstractKatanaSeparator, i, lb, ub, f_tol) = error("Not implemented: Katana.isconstrsat")

"""
    precompute!(sep::AbstractKatanaSeparator, xstar)

Implement this method for a subtype of AbstractKatanaSeparator if your separator might only need to evaluate
certain information once for all constraints using a solution from the internal LP model.

`xstar` is the solution vector from the `KatanaNonlinearModel`'s LP model
"""
precompute!(sep::AbstractKatanaSeparator, xstar) = nothing

### --- First Order Separator Implementation --- #

"""
An implementation of `AbstractKatanaSeparator` for any first-order cutting algorithm.
"""
type KatanaFirstOrderSeparator <: AbstractKatanaSeparator
    sp_cols :: Vector{Vector{Int}} # column indices in Jacobian by row
    sp_col_inds :: Vector{Vector{Int}} # original indices in sparse Jacobian vector

    linear_model :: JuMP.Model # JuMP LP model used by KatanaNonlinearModel
    oracle :: MathProgBase.AbstractNLPEvaluator # evaluator queried as oracle

    num_var :: Int
    num_constr :: Int

    f_tol :: Float64

    # first-order derivative information:
    xstar :: Vector{Float64}
    g     :: Vector{Float64} # vector of constraint values returned by MathProgBase.eval_g()
    jac   :: Vector{Float64} # jacobian sparse matrix

    algo # cutting-plane generator that takes: separator, point and constraint index and returns an AffExpr

    KatanaFirstOrderSeparator(algo) = (s = new(); s.algo = algo; s)
    KatanaFirstOrderSeparator() = KatanaFirstOrderSeparator(linear_oa_cut)
end

# Implements initialize! for KatanaFirstOrderSeparators. This initialises the oracle with necessary features
#  and processes the sparsity structure of the Jacobian to allow faster lookups of non-zero entries by row.
function initialize!(sep          :: KatanaFirstOrderSeparator,
                     linear_model :: JuMP.Model,
                     num_var      :: Int,
                     num_constr   :: Int,
                     f_tol        :: Float64,
                     oracle       :: MathProgBase.AbstractNLPEvaluator)
    sep.linear_model = linear_model

    MathProgBase.initialize(oracle, [:Grad, :Jac, :ExprGraph])
    sep.oracle = oracle

    # map info of Jacobian J_g to allow faster lookups of nonzero entries by row:
    sp_rows, sp_cols = MathProgBase.jac_structure(oracle)
    sep.sp_cols = Vector{Int}[ [] for i=1:num_constr] # map a row to vector of nonzero columns' indices
    sep.sp_col_inds = Vector{Int}[ [] for i=1:num_constr] # map row to vector of indices in sparse matrix
    N = length(sp_rows) # number of sparse entries
    for ind in 1:N
        i,j = sp_rows[ind],sp_cols[ind]
        push!(sep.sp_cols[i], j) # column j is nonzero in row i
        push!(sep.sp_col_inds[i], ind)
    end

    sep.g = zeros(num_constr)
    sep.jac = zeros(N)

    sep.num_var = num_var
    sep.num_constr = num_constr
    sep.f_tol = f_tol
end

# Implements precompute! for KatanaFirstOrderSeparators. This evaluates every constraint at the point xstar and computes
#  and stores the sparse Jacobian at xstar
function precompute!(sep::KatanaFirstOrderSeparator, xstar)
    MathProgBase.eval_jac_g(sep.oracle, sep.jac, xstar) # hopefully variable ordering is consistent with MPB
    MathProgBase.eval_g(sep.oracle, sep.g, xstar) # evaluate constraints

    sep.xstar = xstar # ensure that the x* point matches the gradient information that is precomputed
end

gencuts(sep::KatanaFirstOrderSeparator, xstar, i) = [sep.algo(sep,xstar,i)]

isconstrsat(sep::KatanaFirstOrderSeparator, i, lb, ub) = (sep.g[i] >= lb - sep.f_tol) && (sep.g[i] <= ub + sep.f_tol)

# --- Projection-based Separator Implementation --- #

"""
An implementation of `AbstractKatanaSeparator` that projects a point x_0 to the convex surface and generates a cut
tangent to the constraint.
"""
type KatanaProjectionSeparator <: AbstractKatanaSeparator
    fsep :: KatanaFirstOrderSeparator
    linseps :: Vector{KatanaFirstOrderSeparator}
    projnlps :: Vector{JuMP.Model}
    solver :: MathProgBase.AbstractMathProgSolver
    algo

    num_var :: Int
    eps :: Float64 # epsilon

    xstar :: Vector{Float64}

    function KatanaProjectionSeparator(solver,algo,eps)
        s = new()
        s.algo = algo
        s.solver = solver
        s.eps = eps
        s
    end
    KatanaProjectionSeparator(solver) = KatanaProjectionSeparator(solver,nlp_proj_cut,1e-6)
end

function initialize!(sep          :: KatanaProjectionSeparator,
                     linear_model :: JuMP.Model,
                     num_var      :: Int,
                     num_constr   :: Int,
                     f_tol        :: Float64,
                     oracle       :: MathProgBase.AbstractNLPEvaluator)
    MathProgBase.initialize(oracle, [:ExprGraph])
    sep.num_var = num_var
    sep.projnlps = [JuMP.Model(solver=sep.solver) for i=1:num_constr]
    sep.linseps = [KatanaFirstOrderSeparator() for i=1:num_constr] # separator for each constraint
    for (i,m) in enumerate(sep.projnlps)
        lpvars = [Variable(linear_model,i) for i=1:num_var] # copy bounds on problem variables
        @variable(m, x[i=1:num_var], lowerbound=getlowerbound(lpvars[i]), upperbound=getupperbound(lpvars[i]))
        ex = MathProgBase.constr_expr(oracle, i)
        JuMP.initNLP(m)
        nd, values = JuMP.ReverseDiffSparse.expr_to_nodedata(ex.args[2],m.nlpdata.user_operators)
        nlex = JuMP.NonlinearExprData(nd, values)

        # stolen shamelessly from JuMP's implementation of addNLconstraint:
        op = ex.args[1]
        if op == :(==)
            lb = 0.0
            ub = 0.0
        elseif op == :(<=)
            lb = -Inf
            ub = 0.0
        elseif op == :(>=)
            lb = 0.0
            ub = Inf
        else
            error("unexpected comparison operator")
        end
        
        c = JuMP.NonlinearConstraint(nlex, lb, ub)
        push!(m.nlpdata.nlconstr, c)

        # initialise an NLPEvaluator for each constraint
        d = JuMP.NLPEvaluator(m)
        initialize!(sep.linseps[i], linear_model, num_var, num_constr, f_tol, d)
    end
end

function precompute!(sep::KatanaProjectionSeparator, xstar)
    sep.xstar = xstar
    for fs in sep.linseps
        precompute!(fs, xstar)
    end
end

gencuts(sep::KatanaProjectionSeparator, xstar, i) = sep.algo(sep, xstar, i)

isconstrsat(sep::KatanaProjectionSeparator, i, lb, ub) = isconstrsat(sep.linseps[i], 1, lb, ub)

