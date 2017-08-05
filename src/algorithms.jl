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

# curried function with slurped arguments, will it work?
diff_norm(x0...) = (x...) -> sqrt(sum( (x[i] - x0[i])^2 for i=1:length(x0) ))

function nlp_proj_cut(sep::KatanaProjectionSeparator, a, i::Int)
    m = sep.projnlps[i]

    x = [JuMP.Variable(m, i) for i=1:sep.num_var]
    @NLobjective(m, :Min, sqrt(sum( (x[i] - sep.xstar[i])^2 for i=1:sep.num_var )))

    status = solve(m) # solve projection NLP
    @assert status == :Optimal

    xstar = MathProgBase.getsolution(internalmodel(m))
    println(xstar)
    precompute!(sep.fsep, xstar) # this does have to evaluate every constraint
                                 # in the future, could setup an fsep for every constraint
                                 # using an AbstractNLPEvaluator for each projnlp model
    gencut(sep.fsep, xstar, i)
end
