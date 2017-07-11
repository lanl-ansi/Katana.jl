
type SparseCol
    col :: Int # column index in Jacobian
    ind :: Int # original index in sparse Jacobian vector
end

type KatanaFeatures
    VisData :: Bool
end
KatanaFeatures() = KatanaFeatures(false)

"""
docs go here
"""
type KatanaNonlinearModel <: MathProgBase.AbstractNonlinearModel
    lp_solver    :: MathProgBase.AbstractMathProgSolver
    linear_model :: Union{Void,JuMP.Model}
    oracle       :: MathProgBase.AbstractNLPEvaluator

    status       :: Symbol
    objval       :: Float64

    params       :: KatanaModelParams
    features     :: KatanaFeatures

    num_constr   :: Int64
    num_var      :: Int64
    l_constr     :: Vector{Float64}
    u_constr     :: Vector{Float64}
    l_obj        :: Float64
    u_obj        :: Float64
    objislinear  :: Bool

    nlconstr_ixs :: Vector{Int64} # indices of the NL constraints

    # visualisation logging
    linear_cuts  :: Vector{ConstraintRef{Model,LinearConstraint}}
    lp_sols      :: Vector{Vector{Float64}}

    KatanaNonlinearModel() = new()
end

function KatanaNonlinearModel(lps::MathProgBase.AbstractMathProgSolver, feats::Vector{Symbol}, model_params::KatanaModelParams)
    katana = KatanaNonlinearModel() # don't initialise everything yet
    katana.lp_solver = lps
    katana.status = :None
    katana.objval = NaN

    katana.params = model_params

    katana.features = KatanaFeatures()
    for f in feats
        setfield!(katana.features, f, true)
    end

    katana.nlconstr_ixs = Vector{Int}()
    katana.linear_cuts = Vector{ConstraintRef{Model,LinearConstraint}}()
    katana.lp_sols = Vector{Vector{Float64}}()
    return katana
end

function MathProgBase.NonlinearModel(s::KatanaSolver)
    return KatanaNonlinearModel(s.lp_solver, s.features, s.model_params)
end

# sp_row: vector of SparseCol (nonzero columns in the sparse Jacobian) for a given row
# J: sparse Jacobian matrix, indexed by sp_row
# g: evaluated constraint g_i(a)
# a: point around which we construct the cut
function _constructNewtonCut(m::KatanaNonlinearModel, sp_row::Vector{SparseCol}, J::Vector{Float64}, g::Float64, a)
    # construct the affine expression from sparse gradient:
    #  g'(x) = g_i(a) + (x-a) ⋅ ∇g_i(a)
    v = Vector{JuMP.Variable}()
    coefs = Vector{Float64}()
    b = g
    inner_model = m.linear_model
    for spc in sp_row
        # lookup JuMP variable from column index:
        var = JuMP.Variable(inner_model, spc.col)
        push!(v,var)
        partial = J[spc.ind]
        push!(coefs, partial)
        b += -a[spc.col]*partial
    end
    AffExpr(v, coefs, b) # return an affine expression
end

# add a cut to the internal LP
function _addCut(m::KatanaNonlinearModel, cut::AffExpr, lb::Float64, ub::Float64)
    c = cut.constant
    newconstr = LinearConstraint(cut, lb-c, ub-c)
    cref = JuMP.addconstraint(m.linear_model, newconstr) # add this cut to the LP
    m.features.VisData && push!(m.linear_cuts, cref)
end

function _addEpigraphCut(m::KatanaNonlinearModel, f::Float64, pt)
    ∇f = zeros(m.num_var+1)
    MathProgBase.eval_grad_f(m.oracle, ∇f, pt[1:end-1])
    ∇f[end] = -1
    sp_row = SparseCol[ SparseCol(i,i) for i=1:m.num_var+1 ] # this is a dense row
    y = JuMP.Variable(m.linear_model, m.num_var+1)
    cut = _constructNewtonCut(m, sp_row, ∇f, f, pt)
    _addCut(m, cut, m.l_obj, m.u_obj)
end

_isconstrsat(g, lb, ub, tol) = (g >= lb - tol) && (g <= ub + tol)

function MathProgBase.loadproblem!(
    m::KatanaNonlinearModel,
    num_var::Int, num_constr::Int,
    l_var::Vector{Float64}, u_var::Vector{Float64},
    l_constr::Vector{Float64}, u_constr::Vector{Float64},
    sense::Symbol, d::MathProgBase.AbstractNLPEvaluator)

    m.oracle = d # use this AbstractNLPEvaluator as an oracle that we query for ∇f, g, ∇g
    MathProgBase.initialize(m.oracle, [:Grad,:Jac])

    # setup the internal LP model
    m.linear_model = Model(solver=m.lp_solver)

    inner_lpmod = m.linear_model

    # set up LP
    @variable(inner_lpmod, l_var[i] <= x[i=1:num_var] <= u_var[i])
    @variable(inner_lpmod, m.params.aux_lb <= y <= m.params.aux_ub) # add auxiliary variable
    @objective(inner_lpmod, sense, y)

    # initialise other fields of the KatanaNonlinearModel
    m.num_var = num_var # TODO maybe a better way of keeping track of aux var?
    m.num_constr = num_constr # TODO : should this include constraint on aux var?
    m.l_constr = l_constr
    m.u_constr = u_constr
    # determine bounds on epigraph constraint:
    #   for Max, we have t <= f(x) aka f(x) - t >= 0
    #   for Min, we have t >= f(x) aka f(x) - t <= 0
    m.l_obj, m.u_obj = sense == :Max ? (0.0, Inf) : (-Inf, 0.0)



    # copy linear constraints to linear model by constructing Newton cuts
    #  these Newton cuts recreate the original linear constraint
    J = zeros(m.N) # Jacobian of constraint functions
    g = zeros(m.num_constr) # constraint values
    pt = zeros(m.num_var)
    MathProgBase.eval_jac_g(m.oracle, J, pt) # populate Jacobian
    MathProgBase.eval_g(m.oracle, g, pt) # populate g
    for i=1:num_constr
        if MathProgBase.isconstrlinear(m.oracle, i)
            cut = _constructNewtonCut(m, m.sp_by_row[i], J, g[i], pt)
            _addCut(m, cut, l_constr[i], u_constr[i])
        else
            push!(m.nlconstr_ixs, i) # keep track of these NL constraints for later
        end
    end

    # if linear objective, add it as constraint to LP here
    m.objislinear = MathProgBase.isobjlinear(m.oracle)
    pt = zeros(m.num_var+1)
    if m.objislinear
        println("objective is linear")
        f = MathProgBase.eval_f(m.oracle, pt[1:end-1])
        _addEpigraphCut(m, f, pt)
    else
        println("objective is nonlinear")
    end
end

function MathProgBase.optimize!(m::KatanaNonlinearModel)
    # fixpoint algorithm:
    # 1. run LP solver to compute x*
    # 2. for every unsatisfied non-linear constraint (±f_tol):
    #   3. add cut with specified cutting method
    # 4. check convergence (|g(x) - c| <= f_tol for all g) 

    status = :NotSolved

    J = zeros(m.N) # Jacobian of constraint functions
    g = zeros(m.num_constr) # constraint values
    allsat = false
    iter = 0
    while !allsat && iter <= m.params.iter_cap
        iter += 1
        status = solve(m.linear_model)
        if status == :Unbounded
          # run bounding routine
        elseif status != :Optimal break end

        xstar = MathProgBase.getsolution(internalmodel(m.linear_model))
        m.features.VisData && push!(m.lp_sols, xstar)
        MathProgBase.eval_jac_g(m.oracle , J, xstar[1:end-1]) # hopefully variable ordering is consistent with MPB
        MathProgBase.eval_g(m.oracle, g, xstar[1:end-1]) # evaluate constraints
        allsat = true # base case
        for i in m.nlconstr_ixs # iterate only over NL constraints
#            @assert !MathProgBase.isconstrlinear(m.oracle,i)
            sat = _isconstrsat(g[i], m.l_constr[i], m.u_constr[i], m.params.f_tol)
            if !sat # if constraint not satisfied, add Newton cut
                cut = _constructNewtonCut(m, m.sp_by_row[i], J, g[i], xstar)
                _addCut(m, cut, m.l_constr[i], m.u_constr[i])
            end

            allsat &= sat # loop condition: each constraint must be satisfied
        end
        if !m.objislinear # epigraph constraint is examined separately
            f = MathProgBase.eval_f(m.oracle, xstar[1:end-1]) - xstar[end] # assuming y is last variable
            sat = _isconstrsat(f, m.l_obj, m.u_obj, m.params.f_tol)
            if !sat
                _addEpigraphCut(m, f, xstar)
            end
            allsat &= sat
        end
    end

    println("Katana convergence in $iter iterations.")

    assert(status == :Optimal)
    m.status = status
    return m.status
end


function MathProgBase.setwarmstart!(m::KatanaNonlinearModel, x)
    #TODO, not clear what we can do with x, ignore for now
end

MathProgBase.status(m::KatanaNonlinearModel) = m.status
MathProgBase.getobjval(m::KatanaNonlinearModel) = getobjectivevalue(m.linear_model)

# any auxiliary variables will need to be filtered from this at some point
function MathProgBase.getsolution(m::KatanaNonlinearModel)
  x = MathProgBase.getsolution(internalmodel(m.linear_model))
  x[1:end-1] # aux variable should always be last variable
end


