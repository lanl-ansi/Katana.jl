# Katana.jl Documentation

Katana.jl is a MathProgBase solver for Convex NonLinearPrograms (NLPs).  Katana.jl solves NLPs via the [Extended Cutting-Plane (ECP)](http://epubs.siam.org/doi/10.1137/0108053) method, which combines an Linear Programming solver with a cutting-plane generator to solve Convex NLPs.  Katana.jl is well suited for large-scale Convex NLPs where most of the constraints are linear and the nonlinear constraints are sparse.
