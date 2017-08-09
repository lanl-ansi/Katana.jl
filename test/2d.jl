@testset "2D nonlinear problems" begin
    # these tests only contain two variables (x,y) and are easy to visualize

    # a.k.a. 101_01
    @testset "single convex quadratic constraint, both variables non-zero" begin
        m = Model(solver=solver)

        @variable(m, -2 <= x <= 2)
        @variable(m, -2 <= y <= 2)

        @objective(m, Min, -x-y)
        @NLconstraint(m, x^2 + y^2 <= 1.0)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), -2/sqrt(2), atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 1/sqrt(2), atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 1/sqrt(2), atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 101_02
    @testset "single convex quadratic constraint, one variable non-zero" begin
        m = Model(solver=solver)

        @variable(m, -2 <= x <= 2)
        @variable(m, -2 <= y <= 2)

        @objective(m, Min, -x)
        @NLconstraint(m, x^2 + y^2 <= 1.0)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), -1, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 1, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 0, atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 101_03
    @testset "single convex quadratic constraint, one variable non-zero, max objective" begin
        m = Model(solver=solver)

        @variable(m, -2 <= x <= 2)
        @variable(m, -2 <= y <= 2)

        @objective(m, Max, x)
        @NLconstraint(m, x^2 + y^2 <= 1.0)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 1, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 1, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 0, atol=sol_atol, rtol=sol_rtol)
    end


    # a.k.a. 102_01
    @testset "linear and nonlinear constraints, constraint intersection" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)

        @objective(m, Min, -x)
        @NLconstraint(m, x^2 + y^2 <= 1.0)
        @constraint(m, x + y >= 1.2)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), -0.974165743715913, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 0.974165743715913, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 0.2258342542139504, atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 102_02
    @testset "linear and nonlinear constraints, one binding constraint (linear)" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)

        @objective(m, Min, x+y)
        @NLconstraint(m, x^2 + y^2 <= 1.0)
        @constraint(m, x + y >= 1.2)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 1.2, atol=opt_atol, rtol=opt_rtol)
        ### can't test solution point, here are multiple solutions
        #@test isapprox(getvalue(x), 0.974165743715913, atol=sol_atol, rtol=sol_rtol)
        #@test isapprox(getvalue(y), 0.2258342542139504, atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 102_03
    @testset "linear and nonlinear constraints, one binding constraint (nonlienar)" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)

        @objective(m, Max, x+y)
        @NLconstraint(m, x^2 + y^2 <= 1.0)
        @constraint(m, x + y >= 1.2)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 2/sqrt(2), atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 1/sqrt(2), atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 1/sqrt(2), atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 102_04
    @testset "linear and nonlinear constraints, one binding constraint (quadratic objective)" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)

        @NLobjective(m, Min, x^2+y^2)
        @NLconstraint(m, x^2 + y^2 <= 1.0)
        @constraint(m, x + y >= 1.2)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 0.72, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 0.6, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 0.6, atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 102_05
    @testset "linear and nonlinear constraints, no binding constraints (quadratic objective)" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)

        @NLobjective(m, Min, (x-0.65)^2+(y-0.65)^2)
        @NLconstraint(m, x^2 + y^2 <= 1.0)
        @constraint(m, x + y >= 1.2)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 0.0, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 0.65, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 0.65, atol=sol_atol, rtol=sol_rtol)
    end


    # a.k.a. 103_01
    @testset "nonlinear constraint intersection, one binding constraint (inflection point)" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)

        @objective(m, Min, y)
        @NLconstraint(m, x^2 <= y)
        @NLconstraint(m, -x^2 + 1 >= y)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 0.0, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 0.0, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 0.0, atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 103_02
    @testset "nonlinear constraint intersection, one binding constraint (inflection point)" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)

        @objective(m, Min, -y)
        @NLconstraint(m, x^2 <= y)
        @NLconstraint(m, -x^2 + 1 >= y)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), -1.0, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 0.0, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 1.0, atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 103_03
    @testset "nonlinear constraint intersection, one binding constraint (non-inflection point)" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)

        @objective(m, Min, -x-y)
        @NLconstraint(m, x^2 <= y)
        @NLconstraint(m, -x^2 + 1 >= y)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), -5/4, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 2/4, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 3/4, atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 103_04
    @testset "nonlinear constraint intersection, one binding constraint (non-inflection point)" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)

        @objective(m, Min, x+y)
        @NLconstraint(m, x^2 <= y)
        @NLconstraint(m, -x^2 + 1 >= y)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), -1/4, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), -2/4, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y),  1/4, atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 103_05
    @testset "nonlinear constraint intersection, intersection point" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)

        @objective(m, Min, -x)
        @NLconstraint(m, x^2 <= y)
        @NLconstraint(m, -x^2 + 1 >= y)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), -1/sqrt(2), atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 1/sqrt(2), atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 1/2, atol=sol_atol, rtol=sol_rtol)
    end


    # a.k.a. 104_01
    @testset "redundant nonlinear constraint" begin
        # Note feasible set is inside the circle constraint
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)

        @objective(m, Min, -x)
        @NLconstraint(m, x^2 <= y)
        @NLconstraint(m, -x^2 + 1 >= y)
        @NLconstraint(m, x^2 + (y-0.5)^2 <= 1.0)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), -1/sqrt(2), atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 1/sqrt(2), atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 1/2, atol=sol_atol, rtol=sol_rtol)
    end


    # a.k.a. 105_01
# TODO: test fails with Inf on x[1]
#    @testset "e and log expressions, intersection point" begin
#        m = Model(solver=solver)
#
#        @variable(m, x, start=0.1)
#        @variable(m, y)
#
#        @objective(m, Min, -x-y)
#        @NLconstraint(m, e^(x-2.0) - 0.5 <= y)
#        @NLconstraint(m, log(x) + 0.5 >= y)
#
#        status = solve(m)
#
#        @test status == :Optimal
#        @test isapprox(getobjectivevalue(m), -4.176004405036646, atol=opt_atol, rtol=opt_rtol)
#        @test isapprox(getvalue(x), 2.687422019398147, atol=sol_atol, rtol=sol_rtol)
#        @test isapprox(getvalue(y), 1.488582385638499, atol=sol_atol, rtol=sol_rtol)
#    end

    # a.k.a. 105_02
# TODO: test fails "Adding range constraints not supported yet"
#    @testset "e and log expressions, intersection point" begin
#        m = Model(solver=solver)
#
#        @variable(m, x, start=0.1)
#        @variable(m, y)
#
#        @objective(m, Min, x+y)
#        @NLconstraint(m, e^(x-2.0) - 0.5 <= y)
#        @NLconstraint(m, log(x) + 0.5 >= y)
#
#        status = solve(m)
#
#        @test status == :Optimal
#        @test isapprox(getobjectivevalue(m), 0.16878271368156372, atol=opt_atol, rtol=opt_rtol)
#        @test isapprox(getvalue(x),  0.45538805755556067, atol=sol_atol, rtol=sol_rtol)
#        @test isapprox(getvalue(y), -0.28660534387399694, atol=sol_atol, rtol=sol_rtol)
#    end

    # a.k.a. 105_03
# TODO: test fails "Adding range constraints not supported yet"
#    @testset "e and log expressions, one binding constraint" begin
#        m = Model(solver=solver)
#
#        @variable(m, x, start=0.1)
#        @variable(m, y)
#
#        @objective(m, Min, x-y)
#        @NLconstraint(m, e^(x-2.0) - 0.5 <= y)
#        @NLconstraint(m, log(x) + 0.5 >= y)
#
#        status = solve(m)
#
#        @test status == :Optimal
#        @test isapprox(getobjectivevalue(m), 1/2, atol=opt_atol, rtol=opt_rtol)
#        @test isapprox(getvalue(x), 1, atol=sol_atol, rtol=sol_rtol)
#        @test isapprox(getvalue(y), 1/2, atol=sol_atol, rtol=sol_rtol)
#    end

    # a.k.a. 105_04
# TODO: test fails with Inf on x[1]
#    @testset "e and log expressions, one binding constraint" begin
#        m = Model(solver=solver)
#
#        @variable(m, x, start=0.1)
#        @variable(m, y)
#
#        @objective(m, Min, -x+y)
#        @NLconstraint(m, e^(x-2.0) - 0.5 <= y)
#        @NLconstraint(m, log(x) + 0.5 >= y)
#
#        status = solve(m)
#
#        @test status == :Optimal
#        @test isapprox(getobjectivevalue(m), -3/2, atol=opt_atol, rtol=opt_rtol)
#        @test isapprox(getvalue(x), 2, atol=sol_atol, rtol=sol_rtol)
#        @test isapprox(getvalue(y), 1/2, atol=sol_atol, rtol=sol_rtol)
#    end


    # a.k.a. 106_01
# TODO: test fails for accuracy reasons
#    @testset "cos and sin expressions, convex in given domains, intersection point" begin
#        m = Model(solver=solver)
#
#        @variable(m, -3 <= x <= 3)
#        @variable(m, -1 <= y <= 1)
#
#        @objective(m, Min, -x-y)
#        @NLconstraint(m, sin(-x-1.0) + x/2 + 0.5 <= y)
#        @NLconstraint(m, cos(x-0.5)+ x/4 - 0.5 >= y)
#
#        status = solve(m)
#
#        @test status == :Optimal
#        @test isapprox(getobjectivevalue(m), -1.8572155128552428, atol=opt_atol, rtol=opt_rtol)
#        @test isapprox(getvalue(x), 1.369771397576555, atol=sol_atol, rtol=sol_rtol)
#        @test isapprox(getvalue(y), 0.4874441152786876, atol=sol_atol, rtol=sol_rtol)
#    end

    # a.k.a. 106_02
# TODO: test fails with solver returning Infeasible
#    @testset "cos and sin expressions, convex in given domains, intersection point" begin
#        m = Model(solver=solver)
#
#        @variable(m, -3 <= x <= 3)
#        @variable(m, -1 <= y <= 1)
#
#        @objective(m, Min, x+y)
#        @NLconstraint(m, sin(-x-1.0) + x/2 + 0.5 <= y)
#        @NLconstraint(m, cos(x-0.5)+ x/4 - 0.5 >= y)
#
#        status = solve(m)
#
#        @test status == :Optimal
#        @test isapprox(getobjectivevalue(m), -0.7868226265935826, atol=opt_atol, rtol=opt_rtol)
#        @test isapprox(getvalue(x), -0.5955231764562057, atol=sol_atol, rtol=sol_rtol)
#        @test isapprox(getvalue(y), -0.1912994501373769, atol=sol_atol, rtol=sol_rtol)
#    end


    # a.k.a. 107_01
    @testset "nonlinear objective, nonlinear constraint, no binding constraints" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)

        @NLobjective(m, Min, (x-0.5)^2 + (y-0.5)^2)
        @NLconstraint(m, x^2 + y^2 <= 1)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 0.0, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 0.5, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 0.5, atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 107_02
    @testset "nonlinear objective, nonlinear constraint, binding constraint" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)

        @NLobjective(m, Min, (x-1.0)^2 + (y-1.0)^2)
        @NLconstraint(m, x^2 + y^2 <= 1)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 0.17157287363083387, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 1/sqrt(2), atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 1/sqrt(2), atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 107_03
    @testset "nonlinear objective, nonlinear constraint, binding constraint, different starting point" begin
        m = Model(solver=solver)

        @variable(m, x, start = 1.5)
        @variable(m, y, start = 0.5)

        @NLobjective(m, Min, (x-1.0)^2 + (y-1.0)^2)
        @NLconstraint(m, x^2 + y^2 <= 1)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 0.17157287363083387, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 1/sqrt(2), atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 1/sqrt(2), atol=sol_atol, rtol=sol_rtol)
    end


    # a.k.a. 108_01
    @testset "nonlinear objective,  nonlinear constraint intersection, convex in given domains, no binding constraints" begin
        m = Model(solver=solver)

        @variable(m, x >= 0)
        @variable(m, y >= 0)

        @NLobjective(m, Min, (x-1.0)^2 + (y-0.75)^2)
        @NLconstraint(m, 2*x^2 - 4x*y - 4*x + 4 <= y)
        @NLconstraint(m, y^2 <= -x+2)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 0.0, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 1.00, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 0.75, atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 108_02
# TODO: test fails for accuracy reasons
#    @testset "nonlinear objective,  nonlinear constraint intersection, convex in given domains, intersection point" begin
#        m = Model(solver=solver)
#
#        @variable(m, x >= 0)
#        @variable(m, y >= 0)
#
#        @objective(m, Min, (x-3.0)^2 + y^2)
#        @NLconstraint(m, 2*x^2 - 4x*y - 4*x + 4 <= y)
#        @constraint(m, y^2 <= -x+2)
#
#        status = solve(m)
#
#        @test status == :Optimal
#        @test isapprox(getobjectivevalue(m), 1.5240966871955863, atol=opt_atol, rtol=opt_rtol)
#        @test isapprox(getvalue(x), 1.8344380292075626, atol=sol_atol, rtol=sol_rtol)
#        @test isapprox(getvalue(y), 0.40689308108892147, atol=sol_atol, rtol=sol_rtol)
#    end

    # a.k.a. 108_03
# TODO: test fails for accuracy reasons
#    @testset "nonlinear objective,  nonlinear constraint intersection, convex in given domains, intersection point" begin
#        m = Model(solver=solver)
#
#        @variable(m, x >= 0)
#        @variable(m, y >= 0)
#
#        @objective(m, Min, x^2 + (y-2)^2)
#        @NLconstraint(m, 2*x^2 - 4x*y - 4*x + 4 <= y)
#        @constraint(m, y^2 <= -x+2)
#
#        status = solve(m)
#
#        @test status == :Optimal
#        @test isapprox(getobjectivevalue(m), 0.5927195187027438, atol=opt_atol, rtol=opt_rtol)
#        @test isapprox(getvalue(x), 0.31567986647277146, atol=sol_atol, rtol=sol_rtol)
#        @test isapprox(getvalue(y), 1.2978135998137839, atol=sol_atol, rtol=sol_rtol)
#    end

    # a.k.a. 108_04
# TODO: test fails for accuracy reasons
#    @testset "nonlinear objective,  nonlinear constraint intersection, convex in given domains, one binding constraint" begin
#        m = Model(solver=solver)
#
#        @variable(m, x >= 0)
#        @variable(m, y >= 0)
#
#        @objective(m, Min, x^2 + y^2)
#        @NLconstraint(m, 2*x^2 - 4x*y - 4*x + 4 <= y)
#        @constraint(m, y^2 <= -x+2)
#
#        status = solve(m)
#
#        @test status == :Optimal
#        @test isapprox(getobjectivevalue(m), 0.8112507770394088, atol=opt_atol, rtol=opt_rtol)
#        @test isapprox(getvalue(x), 0.6557120892286371, atol=sol_atol, rtol=sol_rtol)
#        @test isapprox(getvalue(y), 0.6174888121082234, atol=sol_atol, rtol=sol_rtol)
#    end


    # a.k.a. 109_01
# TODO: test fails with Inf on x[1]
#    @testset "nonlinear logarithmic objective, inflection point" begin
#        m = Model(solver=solver)
#
#        @variable(m, x, start=0.1)
#        @variable(m, y, start=0.1)
#
#        @NLobjective(m, Max, log(x))
#        @constraint(m, (y-2)^2 <= -x+2)
#
#        status = solve(m)
#
#        @test status == :Optimal
#        @test isapprox(getobjectivevalue(m), log(2), atol=opt_atol, rtol=opt_rtol)
#        @test isapprox(getvalue(x), 2, atol=sol_atol, rtol=sol_rtol)
#        @test isapprox(getvalue(y), 2, atol=sol_atol, rtol=sol_rtol)
#    end

    # a.k.a. 109_02
# TODO: test fails with Inf on x[1]
#    @testset "nonlinear logarithmic objective, non-inflection point" begin
#        m = Model(solver=solver)
#
#        @variable(m, x, start=0.1)
#        @variable(m, y, start=0.1)
#
#        @NLobjective(m, Max, log(x) + log(y))
#        @constraint(m, (y-2)^2 <= -x+2)
#
#        status = solve(m)
#
#        @test status == :Optimal
#        @test isapprox(getobjectivevalue(m), 1.4853479762665618, atol=opt_atol, rtol=opt_rtol)
#        @test isapprox(getvalue(x), 1.8499011869994715, atol=sol_atol, rtol=sol_rtol)
#        @test isapprox(getvalue(y), 2.387425887570236, atol=sol_atol, rtol=sol_rtol)
#    end

    # a.k.a. 109_03
# TODO: test fails with Inf on x[1]
#    @testset "nonlinear logarithmic objective, non-inflection point" begin
#        m = Model(solver=solver)
#
#        @variable(m, x, start=0.1)
#        @variable(m, y, start=0.1)
#
#        @NLobjective(m, Max, log(x+y))
#        @constraint(m, (y-2)^2 <= -x+2)
#
#        status = solve(m)
#
#        @test status == :Optimal
#        @test isapprox(getobjectivevalue(m), log(7/4 + 5/2), atol=opt_atol, rtol=opt_rtol)
#        @test isapprox(getvalue(x), 7/4, atol=sol_atol, rtol=sol_rtol)
#        @test isapprox(getvalue(y), 5/2, atol=sol_atol, rtol=sol_rtol)
#    end


    # a.k.a. 110_01
    @testset "nonlinear e objective, binding nonlinear constraint, one variable non-zero" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)

        @NLobjective(m, Min, e^(x))
        @NLconstraint(m, x^2 + y^2 <= 1.0)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), e^(-1), atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), -1.0, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y),  0.0, atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 110_02
    @testset "nonlinear e objective, binding nonlinear constraint, both variables non-zero" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)

        @NLobjective(m, Min, e^(x) + e^(y))
        @NLconstraint(m, x^2 + y^2 <= 1.0)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 2*e^(-1/sqrt(2)), atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), -1/sqrt(2), atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), -1/sqrt(2), atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 110_03
    @testset "nonlinear e objective, binding nonlinear constraint, both variables non-zero" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)

        @NLobjective(m, Min, e^(x+y))
        @NLconstraint(m, x^2 + y^2 <= 1.0)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), e^(-2/sqrt(2)), atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), -1/sqrt(2), atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), -1/sqrt(2), atol=sol_atol, rtol=sol_rtol)
    end

end
