@testset "nD nonlinear problems" begin

    # a.k.a. 501_01
    @testset "nD sphere, all variables non-zero" begin
        function nd_shpere(n=2)
            m = Model(solver=solver)

            @variable(m, vars[1:n])

            @objective(m, Min, sum(-x for x in vars))
            @NLconstraint(m, sum(x^2 for x in vars) <= 1.0)

            status = solve(m)

            #println(getobjectivevalue(m))
            #println(getvalue(vars))

            @test status == :Optimal
            #println("$(getobjectivevalue(m)) - $(-n/sqrt(n))")
            @test isapprox(getobjectivevalue(m), -n/sqrt(n), atol=opt_atol, rtol=opt_rtol)
            for x in vars
                #println("$(n) - $(getvalue(x)) - $(1/sqrt(n))")
                @test isapprox(getvalue(x), 1/sqrt(n), atol=sol_atol, rtol=sol_rtol)
            end
        end

        for n in 1:20
            nd_shpere(n)
        end
    end

    # a.k.a. 501_02 (should converge faster than 501_01)
    @testset "nD sphere, all variables non-zero, norm form" begin
        function nd_shpere(n=2)
            m = Model(solver=solver)

            @variable(m, vars[1:n])

            @objective(m, Min, sum(-x for x in vars))
            @NLconstraint(m, sqrt(sum(x^2 for x in vars)) <= 1.0)

            status = solve(m)

            #println(getobjectivevalue(m))
            #println(getvalue(vars))

            @test status == :Optimal
            @test isapprox(getobjectivevalue(m), -n/sqrt(n), atol=opt_atol, rtol=opt_rtol)
            for x in vars
                @test isapprox(getvalue(x), 1/sqrt(n), atol=sol_atol, rtol=sol_rtol)
            end
        end

        for n in 1:20
            nd_shpere(n)
        end
    end

end
