# Katana.jl Documentation

Katana.jl is a MathProgBase solver for Convex NonLinearPrograms (NLPs).  Katana.jl solves NLPs via the [Extended Cutting-Plane (ECP)](http://epubs.siam.org/doi/10.1137/0108053) method, which combines an Linear Programming solver with a cutting-plane generator to solve Convex NLPs.  Katana.jl is well suited for large-scale Convex NLPs where most of the constraints are linear and the nonlinear constraints are sparse.

### Example use

Katana can be used as a solver within a JuMP model. Consider the following non-linear program:

```julia
using Katana, JuMP, GLPKMathProgInterface

# use Katana with default parameters and GLPK as internal LP solver
katana = KatanaSolver(GLPKSolverLP())

m = Model(solver=katana)

# square-root cone constraint, non-linear constraint intersection example:
@variable(m, x, start=0.1)
@variable(m, y, start=0.1)
@variable(m, z)

@objective(m, Min, x+y)
@NLconstraint(m, sqrt(x^2 + y^2) <= z-0.25)
@constraint(m, x^2 + y^2 <= -z + 1)

solve(m)
```

For details on solver parameters, see the [Library documentation](https://lanl-ansi.github.io/Katana.jl/latest/library.html).

### Customising Katana for your use case

Katana is designed with modularity in mind. To that end, although the [default cutting plane algorithm](https://github.com/lanl-ansi/Katana.jl/blob/master/src/algorithms.jl) creates separating hyperplanes by performing a single iteration of [Newton-Raphson](https://en.wikipedia.org/wiki/Newton%27s_method) around the optimal solution found by the linearised model, other separation oracles can be substituted through Katana's [Separators API](https://github.com/lanl-ansi/Katana.jl/blob/master/src/separators.jl).

In addition, you may see benefits to substituting a proprietary linear solver (e.g. Gurobi) for improved stability and speed.

Katana appears to also be well-suited for computing decently tight lower bounds on non-linear objective functions. By specifying the `obj_eps` parameter to the `KatanaSolver`, you can control the solver's stopping criterion to be a 'good enough' objective value.

