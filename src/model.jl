type KatanaFeatures
    VisData :: Bool
end
KatanaFeatures() = KatanaFeatures(false)

"""
The MathProgBase non-linear model for Katana.
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
    if !isempty(filter(x -> !isfinite(x), cut.coeffs))
        Base.warn("Nonlinear constraint or objective likely undefined within domain")
        m.status = :Error
        return
    end
    c = cut.constant
    newconstr = LinearConstraint(cut, lb-c, ub-c)
    cref = JuMP.addconstraint(m.linear_model, newconstr) # add this cut to the LP
    m.numcuts += 1
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
    initialize!(fsep, m.linear_model, num_var+1, num_constr+1, epi_d)
    pt = zeros(num_var+1)
    precompute!(fsep, pt)
    for i=1:num_constr
        if MathProgBase.isconstrlinear(d, i)
            cut = gencut(fsep, pt, (l_constr[i], u_constr[i]), i)
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
        cut = gencut(fsep, pt, (0,0), num_constr+1) # get AffExpr for linear objective
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
            precompute!(fsep, vertex)
            cut = gencut(fsep, vertex, (l_obj, u_obj), m.num_constr)
            round_coefs(cut, m.params.cut_coef_rng)
            _addcut(m, cut, l_obj, u_obj)
        end

        d = EpigraphNLPEvaluator(d, num_var+1, num_constr+1) # wrap d
    end
    m.num_nlconstr = length(m.nlconstr_ixs)

    # initialise the AbstractKatanaSeparator passed to the model
    sep = m.params.separator
    initialize!(sep, m.linear_model, m.num_var, m.num_constr, d)
end

function boundroutine(m::KatanaNonlinearModel, ray)
    # bounding: binary search on unbounded ray until feasibility violated,
    #           at that point, we are outside the constraint surface - make a cut
    #           on the violated constraint(s)
    # NB: this method only works if the feasible region intersects with a sphere of radius 4
    #     around the origin.
    for n=2:1023
        x = (2.0^n)*ray
        allsat = true
        precompute!(m.params.separator, x)
        for i in m.nlconstr_ixs
            sat = isconstrsat(m.params.separator, i, m.l_constr[i], m.u_constr[i], m.params.f_tol)
            if !sat
                cut = gencut(m.params.separator, x, (m.l_constr[i], m.u_constr[i]), i)
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

function print_header()
    @printf("%-10s %-15s %-15s %-20s %-20s %-15s\n", "Iteration", "Total cuts", "Cuts added", "Max constr. viol.", "Avg constr. viol.", "Current cuts")
end

function print_stats(m::KatanaNonlinearModel, iter_lastprnt, cuts_lastprnt, max_viol)
    avg = cuts_lastprnt/(iter_lastprnt*m.num_nlconstr)
    # TODO for now, number of active cuts is same as number of cuts since we don't dynamically remove cuts yet
    @printf("%-10d %-15d %-15d %-20d %-20.2f %-15d\n", m.iter, m.numcuts, cuts_lastprnt, max_viol, avg, m.numcuts)
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
    while status == :Unbounded && i < m.num_var
        ray = MathProgBase.getunboundedray(mpb_lp)
        m.params.log_level > 0 && println("Unbounded ray along: $ray")
        boundroutine(m, ray) # bound along unbounded ray of LP
        m.status == :Error && return m.status
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
    max_viol = 0
    obj_prev = Inf
    while !allsat && m.iter < m.params.iter_cap
        m.iter += 1
        status = solve(m.linear_model, suppress_warnings=true)

        if status != :Optimal # give up
            return m.status = status
        end

        xstar = MathProgBase.getsolution(mpb_lp)

        m.features.VisData && push!(m.lp_sols, xstar)
        precompute!(m.params.separator, xstar)
        allsat = true # base case

        cuts_viol = 0
        for i in m.nlconstr_ixs # iterate only over NL constraints, possibly including epigraph constraint
            sat = isconstrsat(m.params.separator, i, m.l_constr[i], m.u_constr[i], m.params.f_tol)
            if !sat # if constraint not satisfied, call separator API to generate the cut
                cut = gencut(m.params.separator, xstar, (m.l_constr[i], m.u_constr[i]), i)
                round_coefs(cut, m.params.cut_coef_rng)
                _addcut(m, cut, m.l_constr[i], m.u_constr[i])
                m.status == :Error && return m.status
                cuts_viol += 1
            end

            allsat &= sat # loop condition: each constraint must be satisfied
        end
        max_viol = max(max_viol, cuts_viol)
        cuts_lastprnt += cuts_viol

        obj = MathProgBase.getobjval(mpb_lp)
        obj_delta = abs((obj_prev - obj)/obj)
        obj_prev = obj

        if m.params.log_level > 0
            r = m.iter % m.params.log_level
            if r == 0
                (m.iter % (m.params.log_level*50) == 0) && print_header()
                print_stats(m, m.params.log_level, cuts_lastprnt, max_viol)
                cuts_lastprnt = 0
                max_viol = 0
            elseif allsat # print on last iteration also
                print_stats(m, r, cuts_lastprnt, max_viol)
            elseif obj_delta <= m.params.obj_eps
                print_stats(m, m.iter, cuts_lastprnt, max_viol)
            end
        end

        # check convergence of objective value (if obj_eps is <0 ignore this feature)
        if obj_delta <= m.params.obj_eps
            break
        end
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

Returns the number of linear cuts added to the model. Includes linear constraints initially present.
"""
numcuts(m::KatanaNonlinearModel) = m.numcuts

MathProgBase.setwarmstart!(m::KatanaNonlinearModel, x) = fill(0.0, length(x))

MathProgBase.status(m::KatanaNonlinearModel) = m.status
MathProgBase.getobjval(m::KatanaNonlinearModel) = getobjectivevalue(m.linear_model)

# any auxiliary variables will need to be filtered from this at some point
MathProgBase.getsolution(m::KatanaNonlinearModel) = MathProgBase.getsolution(internalmodel(m.linear_model))

MathProgBase.getsolvetime(m::KatanaNonlinearModel) = m.soltime
