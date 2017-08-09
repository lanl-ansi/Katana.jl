@testset "linear constraints only" begin
# these tests are relatively easy
# cut generation is only required if the objective is nonlinear 
# and is not supported natively by the given solver linear constraint solver

    # a.k.a. 001_01
    @testset "linear constraints, closed set" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)

        @objective(m, Min, x)
        @constraint(m, x+y <= 5)
        @constraint(m, 2*x-y <= 3)
        @constraint(m, 3*x+9*y >= -10)
        @constraint(m, 10*x-y >= -20)
        @constraint(m, -x+2*y <= 8)
        @NLconstraint(m, x >= -1000)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), -2.0430107680954848, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), -2.0430107680954848, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), -0.4301075068564087, atol=sol_atol, rtol=sol_rtol)
        #TODO add check that Katana generates no cuts
    end

    # a.k.a. 001_02
    @testset "linear constraints, closed set, quadratic objective, no binding constraints" begin
        m = Model(solver=solver)

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
        @test isapprox(getobjectivevalue(m), 0.0, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 1.0, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 2.0, atol=sol_atol, rtol=sol_rtol)
        #TODO add check that Katana generates no cuts, is solver natively supports the objective
    end

    # a.k.a. 002_01
    @testset "linear constraints, open set" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)

        @objective(m, Min, x+y)
        @constraint(m, 1*x-3*y <= 3)
        @constraint(m, 1*x-5*y <= 0)
        @constraint(m, 3*x+5*y >= 15)
        @constraint(m, 7*x+2*y >= 20)
        @constraint(m, 9*x+1*y >= 20)
        @constraint(m, 3*x+7*y >= 17)
        @NLconstraint(m, x >= -1000)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 3.9655172067026196, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 2.4137930845761546, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 1.5517241221264648, atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 002_02
    @testset "linear constraints, open set, quadratic objective, no binding constraints" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)

        @NLobjective(m, Min, (x-3)^2+(y-2)^2)
        @constraint(m, 1*x-3*y <= 3)
        @constraint(m, 1*x-5*y <= 0)
        @constraint(m, 3*x+5*y >= 15)
        @constraint(m, 7*x+2*y >= 20)
        @constraint(m, 9*x+1*y >= 20)
        @constraint(m, 3*x+7*y >= 17)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 0.0, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 3.0, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 2.0, atol=sol_atol, rtol=sol_rtol)
        #TODO add check that Katana generates no cuts, is solver natively supports the objective
    end
end
