export KatanaSolver

type KatanaSolver <: MathProgBase.AbstractMathProgSolver
    lp_solver::MathProgBase.AbstractMathProgSolver
end

#function KatanaSolver(lp_solver::MathProgBase.AbstractMathProgSolver)
#    return KatanaSolver(lp_solver)
#end
