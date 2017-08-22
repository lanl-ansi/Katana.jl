@testset "3D nonlinear problems" begin
    # these tests only contain two variables (x,y,z) and can be visualized

    # a.k.a. 201_01
    @testset "single convex quadratic constraint, all variables non-zero" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)
        @variable(m, z)

        @objective(m, Min, -(x+y+z))
        @constraint(m, x^2 + y^2 + z^2 <= 1.0)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), -3/sqrt(3), atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 1/sqrt(3), atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 1/sqrt(3), atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(z), 1/sqrt(3), atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 201_02
    @testset "single convex quadratic constraint, one variable non-zero" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)
        @variable(m, z)

        @objective(m, Min, -x)
        @constraint(m, x^2 + y^2 + z^2 <= 1.0)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), -1, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 1, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 0, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(z), 0, atol=sol_atol, rtol=sol_rtol)
    end


    # a.k.a. 202_01
    @testset "intersection nonlinear quadratic constraints, one binding constraint, inflection point" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)
        @variable(m, z)

        @objective(m, Min, -z)
        @constraint(m, x^2 + y^2 <= z)
        @constraint(m, x^2 + y^2 <= -z+1)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), -1, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 0.0, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 0.0, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(z), 1.0, atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 202_02
    @testset "intersection nonlinear quadratic constraints, one binding constraint, inflection point" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)
        @variable(m, z)

        @objective(m, Min, z)
        @constraint(m, x^2 + y^2 <= z)
        @constraint(m, x^2 + y^2 <= -z+1)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 0.0, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 0.0, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 0.0, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(z), 0.0, atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 202_03
    @testset "intersection nonlinear quadratic constraints, one binding constraint, non-inflection point" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)
        @variable(m, z)

        @objective(m, Min, -(x+y+2*z))
        @constraint(m, x^2 + y^2 <= z)
        @constraint(m, x^2 + y^2 <= -z+1)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), -9/4, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 1/4, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 1/4, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(z), 7/8, atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 202_04
    @testset "intersection nonlinear quadratic constraints, one binding constraint, non-inflection point" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)
        @variable(m, z)

        @objective(m, Min, x+y+2*z)
        @constraint(m, x^2 + y^2 <= z)
        @constraint(m, x^2 + y^2 <= -z+1)

        status = solve(m)

        @test status == :Optimal
        # TODO, figure out why ipopt does not ensure 1e-8 on this case
        @test isapprox(getobjectivevalue(m), -1/4, rtol=1e-7)
        @test isapprox(getvalue(x), -1/4, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), -1/4, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(z),  1/8, atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 202_05
    @testset "intersection nonlinear quadratic constraints, one binding constraint, intersection set" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)
        @variable(m, z)

        @objective(m, Min, x+y)
        @constraint(m, x^2 + y^2 <= z)
        @constraint(m, x^2 + y^2 <= -z+1)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), -1.0, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), -1/2, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), -1/2, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(z), 1/2, atol=sol_atol, rtol=sol_rtol)
    end


    # a.k.a. 203_01
    @testset "sqrt cone constraint, nonlinear constraint intersection" begin
        m = Model(solver=solver)

        @variable(m, x, start=0.1)
        @variable(m, y, start=0.1)
        @variable(m, z)

        @objective(m, Min, x+y)
        @NLconstraint(m, sqrt(x^2 + y^2) <= z-0.25)
        @constraint(m, x^2 + y^2 <= -z+1)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), -1/sqrt(2), atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), -sqrt(1/8), atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), -sqrt(1/8), atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(z), 3/4, atol=sol_atol, rtol=sol_rtol)
    end

    # TODO: not supported in JuMP 0.17, try with JuMP 1.0+
    #=
    @testset "sqrt cone constraint (norm form), nonlinear constraint intersection" begin
        m = Model(solver=solver)
        @variable(m, x, start=0.1)
        @variable(m, y, start=0.1)
        @variable(m, z)

        @objective(m, Min, x+y)
        @constraint(m, norm([x,y]) <= z-0.25)
        @constraint(m, x^2 + y^2 <= -z+1)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), -1/sqrt(2), atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), -sqrt(1/8), atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), -sqrt(1/8), atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(z), 3/4, atol=sol_atol, rtol=sol_rtol)
    end
    =#

    # a.k.a. 204_01
    # TODO: test fails with Inf on x[1], issue with initial boundedness of LP
    #=
    @testset "rotated second order cone, nonlinear constraint intersection" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y >= 0)
        @variable(m, z >= 0, start=0.1)

        @objective(m, Min, -y-x)
        @NLconstraint(m, x^2/z <= y)
        @constraint(m, x^2 + y^2 <= -z+1)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), -1.2071067837918394, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 0.353553392657669, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 0.8535533911341705, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(z), 0.14644661317207716, atol=sol_atol, rtol=sol_rtol)
    end
    =#


    # a.k.a. 205_01
    @testset "exponential cone, nonlinear constraint intersection" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y >= 0, start=0.1)
        @variable(m, z)

        @objective(m, Max, y)
        @NLconstraint(m, y*e^(x/y) <= z)
        @NLconstraint(m, y*e^(-x/y) <= z)
        @constraint(m, x^2 + y^2 <= -z+5)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 1.7912878443121907, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 0.0, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 1.7912878443121907, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(z), 1.7912878443121907, atol=sol_atol, rtol=sol_rtol)
    end


    # a.k.a. 206_01
    # TODO: never converges...
    #=
    @testset "power cone, nonlinear constraint intersection" begin
        m = Model(solver=solver)

        # NOTE, starting any of these at 0.0 will segfault libcoinmumps
        @variable(m, x >= 0, start=0.1)
        @variable(m, y >= 0, start=0.1)
        @variable(m, 0 <= z <= 10, start=0.1)

        @objective(m, Max, 2*x+y+z)
        @NLconstraint(m, z <= x^(0.3)*y^(0.7))
        @NLconstraint(m, x <= z^(0.7)*y^(0.3))
        @constraint(m, x^2 + y^2 <= z+1)

        status = solve(m)

        @test status == :Optimal
        # TODO, figure out why ipopt does not ensure 1e-8 on this case
        @test isapprox(getobjectivevalue(m), 4.0, rtol=1e-7)
        @test isapprox(getvalue(x), 1.0, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 1.0, rtol=1e-7)
        @test isapprox(getvalue(z), 1.0, rtol=1e-7)
    end
    =#

    # a.k.a. 210_01
    @testset "nonlinear objective, nonlinear constraint, no binding constraints" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)
        @variable(m, z)

        @objective(m, Min, (x-0.5)^2 + (y-0.5)^2 + (z-0.5)^2)
        @NLconstraint(m, x^2 + y^2 + z^2 <= 1.0)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 0.0, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 1/2, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 1/2, atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(z), 1/2, atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 210_02
    @testset "nonlinear objective, nonlinear constraint, binding constraints" begin
        m = Model(solver=solver)

        @variable(m, x)
        @variable(m, y)
        @variable(m, z)

        @objective(m, Min, (x-1.0)^2 + (y-1.0)^2 + (z-1.0)^2)
        @NLconstraint(m, x^2 + y^2 + z^2 <= 1.0)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 0.535898380052066, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 1/sqrt(3), atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 1/sqrt(3), atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(z), 1/sqrt(3), atol=sol_atol, rtol=sol_rtol)
    end

    # a.k.a. 210_03
    @testset "nonlinear objective, nonlinear constraint, binding constraints, different starting point" begin
        m = Model(solver=solver)

        @variable(m, x, start=1.5)
        @variable(m, y, start=1.0)
        @variable(m, z, start=0.5)

        @objective(m, Min, (x-1.0)^2 + (y-1.0)^2 + (z-1.0)^2)
        @NLconstraint(m, x^2 + y^2 + z^2 <= 1.0)

        status = solve(m)

        @test status == :Optimal
        @test isapprox(getobjectivevalue(m), 0.535898380052066, atol=opt_atol, rtol=opt_rtol)
        @test isapprox(getvalue(x), 1/sqrt(3), atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(y), 1/sqrt(3), atol=sol_atol, rtol=sol_rtol)
        @test isapprox(getvalue(z), 1/sqrt(3), atol=sol_atol, rtol=sol_rtol)
    end

end
