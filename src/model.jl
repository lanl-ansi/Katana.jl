
type KatanaNonlinearModel <: MathProgBase.AbstractNonlinearModel
    lp_solver::MathProgBase.AbstractMathProgSolver
    linear_model::Union{Void,JuMP.Model}
    status::Symbol
    objval::Float64
end

function MathProgBase.NonlinearModel(s::KatanaSolver)
    return KatanaNonlinearModel(s.lp_solver, nothing, :None, NaN)
end


function MathProgBase.loadproblem!(
    m::KatanaNonlinearModel,
    num_var::Int, num_constr::Int,
    l_var::Vector{Float64}, u_var::Vector{Float64},
    l_constr::Vector{Float64}, u_constr::Vector{Float64},
    sense::Symbol, d::MathProgBase.AbstractNLPEvaluator)

    MathProgBase.initialize(d, [:Grad,:Jac,:Hess,:ExprGraph])

    # setup the internal LP model
    m.linear_model = Model(solver=m.lp_solver)

    # add variables
    #@variable(m.linear_model, l_var <= x[1:num_var] <= u_var)

    # by convention, "x" variables can be the original variables and 
    # "y" variables can be auxiliary variables
    # a convention along these lines will help with filtering later

    if MathProgBase.isobjlinear(d)
        # add to model
        println("objective is linear")
    else
        if MathProgBase.isobjquadratic(d) && False # (add check if m.lp_solver can support quadratic obj)
            # add to model
        else
            # add aux variable for objective and add to model
            println("objective is nonlinear")
        end
    end

    for i in 1:num_constr
        if MathProgBase.isconstrlinear(d,i)
            #TODO add to model
            println("constraint $(i) is linear")
            println(MathProgBase.constr_expr(d, i))
            # add constraints
            #@constraint(m.linear_model, MathProgBase.constr_expr(d, i))
        end
    end

end


function MathProgBase.optimize!(m::KatanaNonlinearModel)
    status = solve(m.linear_model)

    if status == :Unbounded
        #run bounding routine
    end

    #TODO add fixpoint algorithm

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
MathProgBase.getsolution(m::KatanaNonlinearModel) = MathProgBase.getsolution(internalmodel(m.linear_model))


