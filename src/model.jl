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
    objislinear  :: Bool

    nlconstr_ixs :: Vector{Int64} # indices of the NL constraints

    # visualisation logging
    linear_cuts  :: Vector{ConstraintRef{Model,LinearConstraint}}
    lp_sols      :: Vector{Vector{Float64}}

    # stats
    iter   :: Int

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

    katana.iter = 0
    return katana
end

function MathProgBase.NonlinearModel(s::KatanaSolver)
    return KatanaNonlinearModel(s.lp_solver, s.features, s.model_params)
end

# add a cut to the internal LP
function _addcut(m::KatanaNonlinearModel, cut::AffExpr, lb::Float64, ub::Float64)
    c = cut.constant
    newconstr = LinearConstraint(cut, lb-c, ub-c)
    cref = JuMP.addconstraint(m.linear_model, newconstr) # add this cut to the LP
    m.features.VisData && push!(m.linear_cuts, cref)
end

function MathProgBase.loadproblem!(
    m::KatanaNonlinearModel,
    num_var::Int, num_constr::Int,
    l_var::Vector{Float64}, u_var::Vector{Float64},
    l_constr::Vector{Float64}, u_constr::Vector{Float64},
    sense::Symbol, d::MathProgBase.AbstractNLPEvaluator)

    # setup the internal LP model
    m.linear_model = Model(solver=m.lp_solver)

    # set up LP
    @variable(m.linear_model, l_var[i] <= x[i=1:num_var] <= u_var[i])
    @variable(m.linear_model, m.params.aux_lb <= y <= m.params.aux_ub) # add auxiliary variable
    @objective(m.linear_model, sense, y)

    # initialise other fields of the KatanaNonlinearModel
    m.num_var = num_var+1 # NB: counts aux variable
    m.num_constr = num_constr+1 # NB: counts epigraph constraint
    m.l_constr = l_constr
    m.u_constr = u_constr
    # determine bounds on epigraph constraint:
    #   for Max, we have t <= f(x) aka f(x) - t >= 0
    #   for Min, we have t >= f(x) aka f(x) - t <= 0
    l_obj, u_obj = sense == :Max ? (0.0, Inf) : (-Inf, 0.0)
    push!(m.l_constr, l_obj)
    push!(m.u_constr, u_obj)

    # initialise the AbstractKatanaSeparator
    sep = m.params.separator
    initialize!(sep, m.linear_model, m.num_var, m.num_constr, d)

    # copy linear constraints to linear model by constructing Newton cuts
    #  these Newton cuts recreate the original linear constraint
    fsep = KatanaFirstOrderSeparator() # with default algo
    initialize!(fsep, m.linear_model, m.num_var, m.num_constr, d)
    pt = zeros(m.num_var)
    precompute!(fsep, pt)
    for i=1:num_constr
        if MathProgBase.isconstrlinear(d, i)
            cut = gencut(fsep, pt, i)
            _addcut(m, cut, l_constr[i], u_constr[i])
        else
            push!(m.nlconstr_ixs, i) # keep track of these NL constraints for later
        end
    end

    # if linear objective, add it as constraint to LP here
    m.objislinear = MathProgBase.isobjlinear(d)
    if m.objislinear
        println("objective is linear")
        cut = gencut(fsep, pt, m.num_constr)
        _addcut(m, cut, l_obj, u_obj)
    else
        println("objective is nonlinear")
        push!(m.nlconstr_ixs, m.num_constr)
    end
end

function MathProgBase.optimize!(m::KatanaNonlinearModel)
    # fixpoint algorithm:
    # 1. run LP solver to compute x*
    # 2. for every unsatisfied non-linear constraint (±f_tol):
    #   3. add cut with specified cutting method
    # 4. check convergence (|g(x) - c| <= f_tol for all g) 

    status = :NotSolved

    allsat = false
    while !allsat && m.iter <= m.params.iter_cap
        m.iter += 1
        status = solve(m.linear_model)
        if status == :Unbounded
          # run bounding routine
        elseif status != :Optimal break end

        xstar = MathProgBase.getsolution(internalmodel(m.linear_model))
        m.features.VisData && push!(m.lp_sols, xstar)
        precompute!(m.params.separator, xstar)
        allsat = true # base case
        for i in m.nlconstr_ixs # iterate only over NL constraints, possibly including epigraph constraint
            sat = isconstrsat(m.params.separator, i, m.l_constr[i], m.u_constr[i], m.params.f_tol)
            if !sat # if constraint not satisfied, call separator API to generate the cut
                # NB: it is the cut-generation algo's responsibility to check if this constraint
                #  is on the objective and handle it accordingly
                cut = gencut(m.params.separator, xstar, i)
                _addcut(m, cut, m.l_constr[i], m.u_constr[i])
            end

            allsat &= sat # loop condition: each constraint must be satisfied
        end
    end

    println("Katana convergence in $(m.iter) iterations.")

    assert(status == :Optimal)
    m.status = status
    return m.status
end

"""
    katana_numiters(m::KatanaNonlinearModel)

Returns the number of iterations taken by the model.
"""
katana_numiters(m::KatanaNonlinearModel) = m.iter

MathProgBase.setwarmstart!(m::KatanaNonlinearModel, x) = fill(NaN, m.num_var-1)

MathProgBase.status(m::KatanaNonlinearModel) = m.status
MathProgBase.getobjval(m::KatanaNonlinearModel) = getobjectivevalue(m.linear_model)
MathProgBase.numconstr(m::KatanaNonlinearModel) = MathProgBase.numconstr(m.linear_model)

# any auxiliary variables will need to be filtered from this at some point
function MathProgBase.getsolution(m::KatanaNonlinearModel)
  x = MathProgBase.getsolution(internalmodel(m.linear_model))
  x[1:end-1] # aux variable should always be last variable
end


