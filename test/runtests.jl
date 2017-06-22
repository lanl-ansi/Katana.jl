using Base.Test

using JuMP

using Ipopt
using GLPKMathProgInterface
using Katana

opt_tol = 1e-8
sol_tol = 1e-8

ipopt = IpoptSolver(print_level=0)

katana = KatanaSolver(GLPKSolverLP())


m = Model()

@variable(m, x)
@variable(m, y)

@objective(m, Min, x)
@constraint(m, x+y <= 5)
@constraint(m, 2*x-y <= 3)
@constraint(m, 3*x+9*y >= -10)
@constraint(m, 10*x-y >= -20)
@constraint(m, -x+2*y <= 8)

# dummy constraint to force NonlinearModel instead of LinearQuadraticModel
@NLconstraint(m, x >= -1000)


setsolver(m, ipopt)
status = solve(m)
@test status == :Optimal

println("Ipopt Solve")
println(getobjectivevalue(m))
println(getvalue(x))
println(getvalue(y))
println("")


setsolver(m, katana)
status = solve(m)
@test status == :Optimal

println("Katana+GLPK Solve")
println(getobjectivevalue(m))
println(getvalue(x))
println(getvalue(y))
println("")

