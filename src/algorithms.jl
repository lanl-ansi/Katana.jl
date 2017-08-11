# this requires a first order separator for which first-order information near a point has been precomputed
#  that means, `a` is not used, but exists for API compatibility
function linear_oa_cuts(sep::KatanaFirstOrderSeparator, a, i::Int)
    # construct the affine expression from sparse gradient:
    #  g_i(a) + (x-a) ⋅ ∇g_i(a)
    v = Vector{JuMP.Variable}()
    coefs = Vector{Float64}()
    b = sep.g[i] # evaluated constraint
    for (col,colix) in zip(sep.sp_cols[i], sep.sp_col_inds[i]) # sparse row
        # lookup JuMP variable from column index:
        var = JuMP.Variable(sep.linear_model, col)
        push!(v,var)
        partial = sep.jac[colix]
        push!(coefs, partial)
        b += -sep.xstar[col]*partial
    end
    [AffExpr(v, coefs, b)] # return an affine expression
end

function _densifygrad(aff::AffExpr, N::Int)
    ∇gi = zeros(N)
    for (v,c) in zip(aff.vars, aff.coeffs)
        ∇gi[v.col] += c # densify gradient
    end
    return ∇gi
end

function orthonormal_oa_cuts(linsep::KatanaFirstOrderSeparator, a, i)
    cuts = Vector{AffExpr}()
    xstar = linsep.xstar
    affs = linear_oa_cuts(linsep, xstar, i)
    ∇g = _densifygrad(affs[1], linsep.num_var)
    U = nullspace(convert(Matrix{Float64}, ∇g')) # orthonormal basis along the tangent hyperplane
    for j=1:size(U,2) # dimension of nullspace is n-1
        u = U[:,j]
        for s=[-1,1] # need to move in both directions
            xprime = xstar + s*u*linsep.eps # move episilon along this vector
            precompute!(linsep, xprime)
            affs = linear_oa_cuts(linsep, xprime, i)
            push!(cuts, affs[1]) # generate a cut at this point
        end
    end
    return cuts
end

function nlp_proj_cut(sep::KatanaProjectionSeparator, a, i::Int)
    m = sep.projnlps[i]
    fsep = sep.linseps[i]

    # linear constraints should be approximated by the FirstOrderSeparator in loadproblem!
    @assert !MathProgBase.isconstrlinear(sep.oracle, i)

    x = [JuMP.Variable(m, i) for i=1:sep.num_var]
    @NLobjective(m, :Min, sum( (x[i] - sep.xstar[i])^2 for i=1:sep.num_var ))

    status = solve(m) # solve projection NLP
    @assert status == :Optimal

    xstar = MathProgBase.getsolution(internalmodel(m))
    println(xstar)
    precompute!(fsep, xstar)

    return gencuts(fsep, xstar, 1)
end
