using Base.Test

using JuMP

using Ipopt
using GLPKMathProgInterface
using Katana

opt_tol = 1e-6
sol_tol = 1e-6

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

ipopt_val = getobjectivevalue(m)
ipopt_x = getvalue(x)
ipopt_y = getvalue(y)

println("Ipopt Solve")
println(ipopt_val)
println(ipopt_x)
println(ipopt_y)
println("")

setsolver(m, katana)
status = solve(m)
@test status == :Optimal

katana_val = getobjectivevalue(m)
katana_x = getvalue(x)
katana_y = getvalue(y)

println("Katana+GLPK Solve")
println(katana_val)
println(katana_x)
println(katana_y)
println("")

# linear program should have identical results
@test abs(katana_val - ipopt_val) <= opt_tol
@test abs(katana_x - ipopt_x) <= sol_tol
@test abs(katana_y - ipopt_y) <= sol_tol

m = Model()

@variable(m, -2 <= x <= 2)
@variable(m, -2 <= y <= 2)

@objective(m, Min, -x-y)
@NLconstraint(m, x^2 + y^2 <= 1.0)

setsolver(m, ipopt)
status = solve(m)
@test status == :Optimal

ipopt_val = getobjectivevalue(m)
ipopt_x = getvalue(x)
ipopt_y = getvalue(y)

println("Ipopt Solve")
println(ipopt_val)
println(ipopt_x)
println(ipopt_y)
println("")

setsolver(m, katana)
status = solve(m)
@test status == :Optimal

katana_val = getobjectivevalue(m)
katana_x = getvalue(x)
katana_y = getvalue(y)

println("Katana+GLPK Solve")
println(katana_val)
println(katana_x)
println(katana_y)
println("")

# cutting-plane approach should produce optimal objective
@test abs(katana_val - ipopt_val) <= opt_tol

