# this requires a first order separator for which first-order information near a point has been precomputed
#  that means, `a` is not used, but exists for API compatibility
function linear_oa_cut(sep::KatanaFirstOrderSeparator, a, i::Int)
    J = sep.jac # sparse jacobian, indexed by sp_row
    g = 0.0 
    sp_row = nothing
    if i > sep.num_constr-1 # add epigraph cut
        g = sep.f # 'constraint' is objective value
        J = sep.∇f # 'sparse Jacobian' will be gradient of f
        J[end] = -1 # aux var gradient always -1
        sp_row = enumerate(1:sep.num_var) # 'sparse row' will just be the entire vector J
    else
        g = sep.g[i] # evaluated constraint
        sp_row = zip(sep.sp_cols[i], sep.sp_col_inds[i]) # sparse row
    end

    # construct the affine expression from sparse gradient:
    #  g'(x) = g_i(a) + (x-a) ⋅ ∇g_i(a)
    v = Vector{JuMP.Variable}()
    coefs = Vector{Float64}()
    b = g
    for (col,colix) in sp_row
        # lookup JuMP variable from column index:
        var = JuMP.Variable(sep.linear_model, col)
        push!(v,var)
        partial = J[colix]
        push!(coefs, partial)
        b += -sep.xstar[col]*partial
    end
    AffExpr(v, coefs, b) # return an affine expression
end


