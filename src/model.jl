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

    # set up LP variables
    @variable(m.linear_model, l_var[i] <= x[i=1:num_var] <= u_var[i])

    # initialise other fields of the KatanaNonlinearModel
    m.num_var = num_var
    m.num_constr = num_constr
    m.l_constr = l_constr
    m.u_constr = u_constr

    # copy linear constraints to linear model by constructing Newton cuts
    #  these Newton cuts recreate the original linear constraint
    # act as if there is an auxiliary variable so that we can generate a 'cut'
    #  on the objective as well, if the objective is linear, in order to recreate
    #  the original linear objective as a constraint
    fsep = KatanaFirstOrderSeparator() # with default algo
    epi_d = EpigraphNLPEvaluator(d, num_var+1, num_constr+1) # wrap d
    initialize!(fsep, m.linear_model, num_var+1, num_constr+1, epi_d)
    pt = zeros(num_var+1)
    precompute!(fsep, pt)
    for i=1:num_constr
        if MathProgBase.isconstrlinear(d, i)
            cut = gencut(fsep, pt, i)
            _addcut(m, cut, l_constr[i], u_constr[i])
        else
            push!(m.nlconstr_ixs, i) # keep track of these NL constraints for later
        end
    end

    # if nonlinear objective, create new variable and bound it by the objective
    m.objislinear = MathProgBase.isobjlinear(d)
    if m.objislinear
        println("objective is linear")
        # this 'cut' approximates the original linear objective exactly
        cut = gencut(fsep, pt, num_constr+1) # get AffExpr for linear objective
        @assert cut.vars[end].col == num_var+1
        pop!(cut.vars) # pop non-existant 'variable'
        pop!(cut.coeffs) # and its coefficient
        JuMP.setobjective(m.linear_model, sense, cut)
    else
        println("objective is nonlinear")

        @variable(m.linear_model, y) # add auxiliary variable
        m.num_var += 1
        @objective(m.linear_model, sense, y)

        # determine bounds on epigraph constraint and add to constraints
        #   for Max, we have t <= f(x) aka f(x) - t >= 0
        #   for Min, we have t >= f(x) aka f(x) - t <= 0
        l_obj, u_obj = sense == :Max ? (0.0, Inf) : (-Inf, 0.0)
        push!(m.l_constr, l_obj)
        push!(m.u_constr, u_obj)
        m.num_constr += 1

        push!(m.nlconstr_ixs, m.num_constr) # add new constraint to list of NL indices

        d = EpigraphNLPEvaluator(d, num_var+1, num_constr+1) # wrap d
    end

    # initialise the AbstractKatanaSeparator passed to the model
    sep = m.params.separator
    initialize!(sep, m.linear_model, m.num_var, m.num_constr, d)
end

function MathProgBase.optimize!(m::KatanaNonlinearModel)
    # fixpoint algorithm:
    # 1. run LP solver to compute x*
    # 2. for every unsatisfied non-linear constraint (±f_tol):
    #   3. add cut with specified cutting method
    # 4. check convergence (|g(x) - c| <= f_tol for all g) 

    status = :NotSolved

    allsat = false
    while !allsat && m.iter < m.params.iter_cap
        m.iter += 1
        status = solve(m.linear_model)
        # TODO if status == :Unbounded run bounding routine
        if status != :Optimal return m.status = status end

        xstar = MathProgBase.getsolution(internalmodel(m.linear_model))
        m.features.VisData && push!(m.lp_sols, xstar)
        precompute!(m.params.separator, xstar)
        allsat = true # base case
        for i in m.nlconstr_ixs # iterate only over NL constraints, possibly including epigraph constraint
            sat = isconstrsat(m.params.separator, i, m.l_constr[i], m.u_constr[i], m.params.f_tol)
            if !sat # if constraint not satisfied, call separator API to generate the cut
                cut = gencut(m.params.separator, xstar, i)
                _addcut(m, cut, m.l_constr[i], m.u_constr[i])
            end

            allsat &= sat # loop condition: each constraint must be satisfied
        end
    end

    println("Katana convergence in $(m.iter) iterations.")

    m.status = status
    return m.status
end

"""
    numiters(m::KatanaNonlinearModel)

Returns the number of iterations taken by the model.
"""
numiters(m::KatanaNonlinearModel) = m.iter

MathProgBase.setwarmstart!(m::KatanaNonlinearModel, x) = fill(0.0, length(x))

MathProgBase.status(m::KatanaNonlinearModel) = m.status
MathProgBase.getobjval(m::KatanaNonlinearModel) = getobjectivevalue(m.linear_model)
MathProgBase.numconstr(m::KatanaNonlinearModel) = MathProgBase.numconstr(m.linear_model)

# any auxiliary variables will need to be filtered from this at some point
MathProgBase.getsolution(m::KatanaNonlinearModel) = MathProgBase.getsolution(internalmodel(m.linear_model))


