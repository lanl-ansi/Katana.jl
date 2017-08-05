
"""
Wrapper for any `MathProgBase.AbstractNLPEvaluator` instance. Acts as an NLP evaluator for an NLP that has been transformed
into epigraph form by treating the NL objective as if it were a constraint.
"""
type EpigraphNLPEvaluator <: MathProgBase.AbstractNLPEvaluator
    nlpeval :: MathProgBase.AbstractNLPEvaluator # wrapped evaluator
    num_var :: Int # number of variables including any auxiliary variable
    num_constr :: Int # number of constraints including any constraint on objective

    ∇f :: Vector{Float64} # don't bother reallocating space for this on each call to eval_jac_g
end

EpigraphNLPEvaluator(d, num_var, num_constr) = EpigraphNLPEvaluator(d,num_var,num_constr,zeros(num_var))

# pass-through
MathProgBase.isobjlinear(d::EpigraphNLPEvaluator) = MathProgBase.isobjlinear(d.nlpeval)
MathProgBase.isobjquadratic(d::EpigraphNLPEvaluator) = MathProgBase.isobjquadratic(d.nlpeval)
MathProgBase.isconstrlinear(d::EpigraphNLPEvaluator) = MathProgBase.isconstrlinear(d.nlpeval)
MathProgBase.obj_expr(d::EpigraphNLPEvaluator) = MathProgBase.obj_expr(d.nlpeval) # TODO technically objective is f(x) = x_n

# TODO implement any feature requested in order to be a true pass-through to the underlying NLP evaluator
MathProgBase.features_available(d::EpigraphNLPEvaluator) = [:Grad, :Jac, :ExprGraph]

function MathProgBase.initialize(d::EpigraphNLPEvaluator, requested_features)
    for feat in requested_features
        if !(feat in MathProgBase.features_available(d) && feat in MathProgBase.features_available(d.nlpeval))
            error("Unsupported feature $feat")
        end
    end
    MathProgBase.initialize(d.nlpeval, requested_features)
end

# TODO for these, technically f is just the auxiliary variable
MathProgBase.eval_f(d::EpigraphNLPEvaluator, x) = MathProgBase.eval_f(d.nlpeval, x[1:end-1]) - x[end]
function MathProgBase.eval_grad_f(d::EpigraphNLPEvaluator, g, x)
    MathProgBase.eval_grad_f(d.nlpeval, g, x[1:end-1])
    g[end] = -1 # df/dt = -1 since we've constraint is f(x) - t
end

# g should have num_constr elements, x should have num_var elements
function MathProgBase.eval_g(d::EpigraphNLPEvaluator, g, x)
    MathProgBase.eval_g(d.nlpeval, g, x[1:end-1])
    g[end] = MathProgBase.eval_f(d,x)
end

# append (non-sparse) row of jacobian for the epigraph constraint
#  if there is one
function MathProgBase.jac_structure(d::EpigraphNLPEvaluator)
    jstruc = MathProgBase.jac_structure(d.nlpeval)
    append!(jstruc[1], fill(d.num_constr, d.num_var)) # add num_var new elements to end
    append!(jstruc[2], collect(1:d.num_var))
    jstruc
end

# if d.has_epigraph:
# J should have the size N+num_var where N is the original number of nonzero elements in the jacobian,
#  before any epigraph constraint was added
function MathProgBase.eval_jac_g(d::EpigraphNLPEvaluator, J, x)
    MathProgBase.eval_jac_g(d.nlpeval, J, x[1:end-1])
    MathProgBase.eval_grad_f(d, d.∇f, x)
    J[end-d.num_var+1:end] = d.∇f # set last num_var elements to gradient of f
end

