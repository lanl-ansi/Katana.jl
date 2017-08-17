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
    num_nlconstr :: Int

    # TODO: update visualisation stuff
    cut_pool     :: Vector{Int} # map constraint index to age (no. of iterations in a row marked inactive)

    # visualisation logging
    linear_cuts  :: Vector{ConstraintRef{Model,LinearConstraint}}
    lp_sols      :: Vector{Vector{Float64}}

    # stats
    iter   :: Int
    numcuts :: Int
    soltime :: Float64

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
    katana.cut_pool= Vector{Int}()
    katana.linear_cuts = Vector{ConstraintRef{Model,LinearConstraint}}()
    katana.lp_sols = Vector{Vector{Float64}}()

    katana.iter = 0
    katana.numcuts = 0
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
    m.numcuts += 1
    push!(m.cut_pool, 0)
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
    vertex = fill(NaN, num_var) # may need a vertex of the variables' bounds polytope
    status = solve(m.linear_model, suppress_warnings=true)
    if status == :Optimal
        vertex = MathProgBase.getsolution(internalmodel(m.linear_model))
    end

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
    initialize!(fsep, m.linear_model, num_var+1, num_constr+1, m.params.f_tol, epi_d)
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
        m.params.log_level > 0 && println("objective is linear")
        # this 'cut' approximates the original linear objective exactly
        cut = gencut(fsep, pt, num_constr+1) # get AffExpr for linear objective
        @assert cut.vars[end].col == num_var+1
        pop!(cut.vars) # pop non-existant 'variable'
        pop!(cut.coeffs) # and its coefficient
        JuMP.setobjective(m.linear_model, sense, cut)
    else
        m.params.log_level > 0 && println("objective is nonlinear")

        @variable(m.linear_model, y) # add auxiliary variable
        m.num_var += 1
        @objective(m.linear_model, sense, y)

        # determine bounds on epigraph constraint and add to constraints
        #   for Max, we have y <= f(x) aka f(x) - y >= 0
        #   for Min, we have y >= f(x) aka f(x) - y <= 0
        l_obj, u_obj = sense == :Max ? (0.0, Inf) : (-Inf, 0.0)
        push!(m.l_constr, l_obj)
        push!(m.u_constr, u_obj)
        m.num_constr += 1
        push!(m.nlconstr_ixs, m.num_constr) # add new constraint to list of NL indices

        # we can actually introduce a "decent" bound on t by exploiting the
        # convexity of f, assuming all other problem variables x are upper- and lower-bounded
        # 1. solve LP consisting of variable bounds on x, no objective, to compute vertex x*
        #       - this was done earlier, before adding linear constraints
        # 2. at x*, compute f(x*) and a hyperplane tangent to that point
        #       - this is of course done by our cutting-plane algorithm!
        if any(isnan, vertex)
            Base.warn("Problem variables insufficiently bounded!")
        else
            push!(vertex, MathProgBase.eval_f(d, vertex))
            cut = gencut(fsep, vertex, m.num_constr)
            round_coefs(cut, m.params.cut_coef_rng)
            _addcut(m, cut, l_obj, u_obj)
        end

        d = EpigraphNLPEvaluator(d, num_var+1, num_constr+1) # wrap d
    end
    m.num_nlconstr = length(m.nlconstr_ixs)

    # initialise the AbstractKatanaSeparator passed to the model
    sep = m.params.separator
    initialize!(sep, m.linear_model, m.num_var, m.num_constr, m.params.f_tol, d)
end

function boundroutine(m::KatanaNonlinearModel, ray)
    # bounding: binary search on unbounded ray until feasibility violated,
    #           at that point, we are outside the constraint surface - make a cut
    #           on the violated constraint(s)
    for n=2:1023
        x = (2.0^n)*ray
        allsat = true
        precompute!(m.params.separator, x)
        for i in m.nlconstr_ixs
            sat = isconstrsat(m.params.separator, i, m.l_constr[i], m.u_constr[i])
            if !sat
                cut = gencut(m.params.separator, (m.l_constr[i], m.u_constr[i]), i)
                round_coefs(cut, m.params.cut_coef_rng)
                _addcut(m, cut, m.l_constr[i], m.u_constr[i])
            end
            allsat &= sat
        end

        if !allsat break end # stop searching in this direction
    end
end

# determine largest coefficient and set all coefs more than cut_coef_rng lower to 0
function round_coefs( cut::AffExpr, cut_coef_rng::Float64)
    max_coef = maximum(cut.coeffs)
    for i=1:length(cut.coeffs)
        if cut.coeffs[i] + cut_coef_rng < max_coef
            cut.coeffs[i] = 0.0
        end
    end
end

# remove cuts that have been inactive for more than a certain number of iterations
# returns the set of inactive constraints on this iteration
function remove_old_cuts(m::KatanaNonlinearModel)
    mpbmod = internalmodel(m.linear_model)

    # enumerate inactive cuts on this iteration
    dualsol = MathProgBase.getconstrduals(mpbmod)
    inactive_constrs = Set{Int}()
    old_cuts = Vector{Int}()

    for (c,s) in enumerate(dualsol)
        if s == 0
            push!(inactive_constrs, c) # keep track of these
            m.cut_pool[c] += 1 # increment age of this cut
            if m.cut_pool[c] > m.params.max_cut_age
                push!(old_cuts, c)
                m.numcuts -= 1
            end
        else
            m.cut_pool[c] = 0
        end
    end

    # remove these cuts from model and cut pool
    deleteat!(m.linear_model.linconstr, old_cuts)
    deleteat!(m.cut_pool, old_cuts)

    # now we can modify the MPB model
    !isempty(old_cuts) && MathProgBase.delconstrs!(mpbmod, old_cuts)
    inactive_constrs
end

function print_header()
    @printf("%-10s %-15s %-15s %-20s %-20s %-15s\n", "Iteration", "Current cuts", "Cuts added", "Max constr. viol.", "Avg constr. viol.", "Active cuts")
end

function print_stats(m::KatanaNonlinearModel, iter_lastprnt, cuts_lastprnt, max_viol, inactive_cuts)
    avg = cuts_lastprnt/(iter_lastprnt*m.num_nlconstr)
    active_cuts = m.numcuts - inactive_cuts
    # TODO for now, number of active cuts is same as number of cuts since we don't dynamically remove cuts yet
    @printf("%-10d %-15d %-15d %-20d %-20.4f %-15d\n", m.iter, m.numcuts, cuts_lastprnt, max_viol, avg, active_cuts)
end

function MathProgBase.optimize!(m::KatanaNonlinearModel)
    # fixpoint algorithm:
    # 1. run LP solver to compute x*
    # 2. for every unsatisfied non-linear constraint (±f_tol):
    #   3. add cut with specified cutting method
    # 4. check convergence (|g(x) - c| <= f_tol for all g) 

    # presolve: resolve initially-unbounded LP
    start = time()
    status = solve(m.linear_model,suppress_warnings=true)
    if status == :Unbounded
        Base.warn("Automatically bounding unbounded LP")
    end

    mpb_lp = internalmodel(m.linear_model)
    i = 0
    while status == :Unbounded && i < max(m.params.presolve_cap, m.num_var)
        ray = MathProgBase.getunboundedray(mpb_lp)
        m.params.log_level > 0 && println("Unbounded ray along: $ray")
        boundroutine(m, ray) # bound along unbounded ray of LP
        status = solve(m.linear_model,suppress_warnings=true)
        i += 1
    end

    if status == :Unbounded # if it's still unbounded
        Base.warn("Katana could not resolve unbounded LP")
        return m.status = status
    end

    if m.params.log_level > 0
        print_header()
    end

    allsat = false
    cuts_lastprnt = 0 # number of cuts since last printout
    total_cuts = m.numcuts
    max_viol = 0
    obj_prev = Inf
    while !allsat && m.iter < m.params.iter_cap
        m.iter += 1
        status = solve(m.linear_model, suppress_warnings=true)
        num_inactive_cuts = length(remove_old_cuts(m))

        if status != :Optimal # give up
            return m.status = status
        end

        xstar = MathProgBase.getsolution(mpb_lp)

        m.features.VisData && push!(m.lp_sols, xstar)
        precompute!(m.params.separator, xstar)
        allsat = true # base case

        cuts_viol = 0
        for i in m.nlconstr_ixs # iterate only over NL constraints, possibly including epigraph constraint
            sat = isconstrsat(m.params.separator, i, m.l_constr[i], m.u_constr[i])
            if !sat # if constraint not satisfied, call separator API to generate the cut
                cut = gencut(m.params.separator, (m.l_constr[i], m.u_constr[i]), i)
                round_coefs(cut, m.params.cut_coef_rng)
                _addcut(m, cut, m.l_constr[i], m.u_constr[i])
                cuts_viol += 1
            end

            allsat &= sat # loop condition: each constraint must be satisfied
        end
        max_viol = max(max_viol, cuts_viol)
        cuts_lastprnt += cuts_viol
        total_cuts += cuts_viol

        if m.params.log_level > 0
            r = m.iter % m.params.log_level
            if r == 0
                (m.iter % (m.params.log_level*50) == 0) && print_header()
                print_stats(m, m.params.log_level, cuts_lastprnt, max_viol, num_inactive_cuts)
                cuts_lastprnt = 0
                max_viol = 0
            elseif allsat # print on last iteration also
                print_stats(m, r, cuts_lastprnt, max_viol, num_inactive_cuts)
            end
        end

        # check convergence of objective value (if obj_eps is <0 ignore this feature)
        obj = MathProgBase.getobjval(mpb_lp)
        obj_delta = abs((obj_prev - obj)/obj)
        if obj_delta <= m.params.obj_eps
            break
        end
        obj_prev = obj
    end

    m.soltime = time() - start

    if m.iter >= m.params.iter_cap
        status = :UserLimit
    end

    m.status = status
    return m.status
end

"""
    numiters(m::KatanaNonlinearModel)

Returns the number of iterations taken by the model.
"""
numiters(m::KatanaNonlinearModel) = m.iter

"""
    numcuts(m::KatanaNonlinearModel)

Returns the number of cuts added to the model
"""
numcuts(m::KatanaNonlinearModel) = m.numcuts

MathProgBase.setwarmstart!(m::KatanaNonlinearModel, x) = fill(0.0, length(x))

MathProgBase.status(m::KatanaNonlinearModel) = m.status
MathProgBase.getobjval(m::KatanaNonlinearModel) = getobjectivevalue(m.linear_model)

# any auxiliary variables will need to be filtered from this at some point
MathProgBase.getsolution(m::KatanaNonlinearModel) = MathProgBase.getsolution(internalmodel(m.linear_model))

MathProgBase.getsolvetime(m::KatanaNonlinearModel) = m.soltime
