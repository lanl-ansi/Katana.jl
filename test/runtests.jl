using Logging
# suppress warnings during testing
Logging.configure(level=ERROR)

using Base.Test
using JuMP

using Ipopt
using GLPKMathProgInterface
if Pkg.installed("Gurobi") == nothing
    using Gurobi
end
using Katana

opt_tol = 1e-2
sol_tol = 1e-2

ipopt = IpoptSolver(print_level=0)

katana = KatanaSolver(GLPKSolverLP(), log_level=0)
#katana = KatanaSolver(GurobiSolver(OutputFlag=0))

solver = katana
# useful for makeing sure the tests are correct
#solver = ipopt 


include("basic.jl")

include("lpqp.jl")

include("2d.jl")

include("3d.jl")

include("misc.jl")

