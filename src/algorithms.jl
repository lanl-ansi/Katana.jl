# this requires a first order separator for which first-order information near a point has been precomputed
#  that means, `a` is not used, but exists for API compatibility
function linear_oa_cut(sep::KatanaFirstOrderSeparator, a, i::Int)
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
    AffExpr(v, coefs, b) # return an affine expression
end

function _densifygrad(aff::AffExpr, N::Int)
    ∇gi = zeros(N)
    for (v,c) in zip(aff.vars, aff.coeffs)
        ∇gi[v.col] += c # densify gradient
    end
    return ∇gi
end


function nlp_proj_cut(sep::KatanaProjectionSeparator, a, i::Int)
    m = sep.projnlps[i]
    fsep = sep.linseps[i]

    x = [JuMP.Variable(m, i) for i=1:sep.num_var]
    @NLobjective(m, :Min, sqrt(sum( (x[i] - sep.xstar[i])^2 for i=1:sep.num_var )))

    status = solve(m) # solve projection NLP
    @assert status == :Optimal

    xstar = MathProgBase.getsolution(internalmodel(m))
    println(xstar)
    precompute!(fsep, xstar) # this does have to evaluate every constraint
                                 # in the future, could setup an fsep for every constraint
                                 # using an AbstractNLPEvaluator for each projnlp model
    ∇g = _densifygrad(gencut(fsep, xstar, 1), sep.num_var)
    U = nullspace(convert(Matrix{Float64}, ∇g')) # orthonormal basis along the tangent hyperplane
    cuts = Vector{AffExpr}()
    for i=1:size(U,2) # dimension of nullspace is n-1
        u = U[:,i]
        for s=[-1,1] # need to move in both directions
            xprime = xstar + s*u*sep.eps # move episilon along this vector
            precompute!(fsep, xprime)
            push!(cuts, gencut(fsep, xprime, 1)) # generate a cut at this point
        end
    end

    return cuts
end
