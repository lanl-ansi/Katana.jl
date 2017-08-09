@testset "basic tests" begin
    # very simple tests, just to perform basic code correctness checks

    @testset "pass through of linear model components" begin
        m = Model()

        @variable(m, x)
        @variable(m, y)

        @objective(m, Min, x)
        @constraint(m, x+y <= 5)
        @constraint(m, 2*x-y <= 3)
        @constraint(m, 3*x+9*y >= -10)
        @constraint(m, 10*x-y >= -20)
        @constraint(m, -x+2*y <= 8)
        @NLconstraint(m, x >= -1000) # force NonlinearModel

        setsolver(m, ipopt)
        status = solve(m)
        @test status == :Optimal

        ipopt_val = getobjectivevalue(m)
        ipopt_x = getvalue(x)
        ipopt_y = getvalue(y)

        #println("Ipopt Solve: $ipopt_val, $ipopt_x, $ipopt_y")
        #println("")

        setsolver(m, katana)
        status = solve(m)
        @test status == :Optimal

        katana_val = getobjectivevalue(m)
        katana_x = getvalue(x)
        katana_y = getvalue(y)

        #println("Katana+GLPK Solve: $katana_val, $katana_x, $katana_y")
        #println("")

        @test isapprox(katana_val, ipopt_val, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(katana_x, ipopt_x, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(katana_y, ipopt_y, atol=sol_atol, rtol=sol_rtol)
    end

    @testset "fixpoint convergence" begin
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

        #println("Ipopt Solve: $ipopt_val, $ipopt_x, $ipopt_y")
        #println("")

        setsolver(m, katana)
        status = solve(m)
        @test status == :Optimal

        katana_val = getobjectivevalue(m)
        katana_x = getvalue(x)
        katana_y = getvalue(y)

        #println("Katana+GLPK Solve: $katana_val, $katana_x, $katana_y")
        #println("")

        @test isapprox(katana_val, ipopt_val, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(katana_x, ipopt_x, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(katana_y, ipopt_y, atol=sol_atol, rtol=sol_rtol)
    end

    @testset "lifting of nonlinear objective function" begin
        m = Model(solver=katana)

        @variable(m, x)
        @variable(m, y)

        @NLobjective(m, Min, (x-1)^2 + (y-2)^2)
        @constraint(m, x+y <= 5)
        @constraint(m, 2*x-y <= 3)
        @constraint(m, 3*x+9*y >= -10)
        @constraint(m, 10*x-y >= -20)
        @constraint(m, -x+2*y <= 8)

        status = solve(m)
        @test status == :Optimal

        #println("Katana+GLPK Solve: $(getobjectivevalue(m)), $(getvalue(x)), $(getvalue(y))")
        #println("")

        @test isapprox(getobjectivevalue(m), 0.0, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 1.0, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 2.0, atol=sol_atol, rtol=sol_rtol)
    end
end
