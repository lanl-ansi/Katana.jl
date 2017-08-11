using Logging
# suppress warnings during testing
Logging.configure(level=ERROR)

using Base.Test
using JuMP

using Ipopt
using GLPKMathProgInterface
if Pkg.installed("Gurobi") != nothing
    using Gurobi
end

using Katana

opt_rtol = 1e-6
opt_atol = 1e-6

sol_rtol = 1e-3
sol_atol = 1e-3

ipopt = IpoptSolver(print_level=0)

#separator = KatanaFirstOrderSeparator()
separator = KatanaProjectionSeparator(ipopt,orthonormal_cuts=true)
katana = KatanaSolver(GLPKSolverLP(), log_level=0, separator=separator)

# useful for debugging algorithm correctness
#katana = KatanaSolver(GurobiSolver(OutputFlag=0))

solver = katana
# useful for makeing sure the tests are correct
#solver = ipopt 


include("basic.jl")

include("lpqp.jl")

include("2d.jl")

include("3d.jl")

include("misc.jl")

