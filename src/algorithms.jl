# this requires a first order separator for which first-order information near a point has been precomputed
#  that means, `a` is not used, but exists for API compatibility
function linear_oa_cut(sep::KatanaFirstOrderSeparator, a, i::Int)
    # construct the affine expression from sparse gradient:
    #  g'(x) = g_i(a) + (x-a) ⋅ ∇g_i(a)
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

