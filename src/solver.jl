export KatanaSolver

"""
docs go here
"""
type KatanaSolver <: MathProgBase.AbstractMathProgSolver
    lp_solver::MathProgBase.AbstractMathProgSolver
end

#function KatanaSolver(lp_solver::MathProgBase.AbstractMathProgSolver)
#    return KatanaSolver(lp_solver)
#end

# this bridge should make lp/qp models act like nlp models
function MathProgBase.LinearQuadraticModel(s::KatanaSolver)
    MathProgBase.NonlinearToLPQPBridge(MathProgBase.NonlinearModel(s))
end
