# Katana.jl
Katana.jl is a MathProgBase solver for Convex NLPs.  The Katana.jl solves NLP using an LP solver and generating cutting-planes.  Katana.jl is suited for large-scale Convex NLPs, where the most of the constraints are linear and the nonlinear constraints are sparse.


## Installation
Install via `Pkg.clone("git@github.com:lanl-ansi/Katana.jl.git")`

Test via `Pkg.test("Katana")`

Build docs by running `julia make.jl` in the docs directory


## License
This code is provided under a BSD license as part of the Polyhedral Approximation in Julia: Automatic Reformulations for InTeger Optimization (PAJARITO) project, LA-CC-15-088.
