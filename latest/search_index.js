var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Katana.jl-Documentation-1",
    "page": "Home",
    "title": "Katana.jl Documentation",
    "category": "section",
    "text": "Katana.jl is a MathProgBase solver for Convex NonLinearPrograms (NLPs).  Katana.jl solves NLPs via the Extended Cutting-Plane (ECP) method, which combines an Linear Programming solver with a cutting-plane generator to solve Convex NLPs.  Katana.jl is well suited for large-scale Convex NLPs where most of the constraints are linear and the nonlinear constraints are sparse."
},

{
    "location": "index.html#Example-use-1",
    "page": "Home",
    "title": "Example use",
    "category": "section",
    "text": "Katana can be used as a solver within a JuMP model. Consider the following non-linear program:using Katana, JuMP, GLPKMathProgInterface\n\n# use Katana with default parameters and GLPK as internal LP solver\nkatana = KatanaSolver(GLPKSolverLP())\n\nm = Model(solver=katana)\n\n# square-root cone constraint, non-linear constraint intersection example:\n@variable(m, x, start=0.1)\n@variable(m, y, start=0.1)\n@variable(m, z)\n\n@objective(m, Min, x+y)\n@NLconstraint(m, sqrt(x^2 + y^2) <= z-0.25)\n@constraint(m, x^2 + y^2 <= -z + 1)\n\nsolve(m)For details on solver parameters, see the Library documentation."
},

{
    "location": "index.html#Customising-Katana-for-your-use-case-1",
    "page": "Home",
    "title": "Customising Katana for your use case",
    "category": "section",
    "text": "Katana is designed with modularity in mind. To that end, although the default cutting plane algorithm creates separating hyperplanes by performing a single iteration of Newton-Raphson around the optimal solution found by the linearised model, other separation oracles can be substituted through Katana's Separators API.In addition, you may see benefits to substituting a proprietary linear solver (e.g. Gurobi) for improved stability and speed.Katana appears to also be well-suited for computing decently tight lower bounds on non-linear objective functions. By specifying the obj_eps parameter to the KatanaSolver, you can control the solver's stopping criterion to be a 'good enough' objective value.Consult the User Manual for more."
},

{
    "location": "manual.html#",
    "page": "Manual",
    "title": "Manual",
    "category": "page",
    "text": ""
},

{
    "location": "manual.html#Katana-User-Manual-1",
    "page": "Manual",
    "title": "Katana User Manual",
    "category": "section",
    "text": "A more in-depth guide on customising, tweaking and extending Katana."
},

{
    "location": "manual.html#Implementing-Custom-Separators-1",
    "page": "Manual",
    "title": "Implementing Custom Separators",
    "category": "section",
    "text": "The Separator API is intended to abstract the logic of generating separating hyperplanes into two parts:What stateful information is required? For example, the KatanaFirstOrderSeparator stores evaluated constraint values, evaluated constraint Jacobians, Jacobian sparsity structure, etc.\nHow is that information used to generate the cutting plane? In the default Newton-based algorithm, the hyperplane is constructed as a first-order expansion around a point.Implementing a custom separator requires a distinction between the _information_ required and the _method_ used. It may not be necessary to create a new subclass of AbstractKatanaSeparator if, for example, the separation oracle only intends to use first-order information. Consider a method of generating cutting planes that involves starting at the LP optimal point x^* and performing gradient descent until within some epsilon of the constraint surface. Ideally, this would require another function that can be used as the algo of a KatanaFirstOrderSeparator.On the other hand, if the custom separator requires second-order constraint values, such as the Hessian, it would be necessary to implement a new subclass of AbstractKatanaSeparator. It is then up to that implementation if the actual algorithm called by gencut is itself modular."
},

{
    "location": "manual.html#Expressing-models-in-JuMP-for-better-performance-with-Katana-1",
    "page": "Manual",
    "title": "Expressing models in JuMP for better performance with Katana",
    "category": "section",
    "text": "TODO."
},

{
    "location": "library.html#",
    "page": "Library",
    "title": "Library",
    "category": "page",
    "text": ""
},

{
    "location": "library.html#Katana.AbstractKatanaSeparator",
    "page": "Library",
    "title": "Katana.AbstractKatanaSeparator",
    "category": "Type",
    "text": "See the User Manual for an explanation of the motivation behind this class.\n\n\n\n"
},

{
    "location": "library.html#Katana.KatanaFirstOrderSeparator",
    "page": "Library",
    "title": "Katana.KatanaFirstOrderSeparator",
    "category": "Type",
    "text": "An implementation of AbstractKatanaSeparator for any first-order cutting algorithm.\n\n\n\n"
},

{
    "location": "library.html#Katana.KatanaSolver",
    "page": "Library",
    "title": "Katana.KatanaSolver",
    "category": "Type",
    "text": "A solver for convex NLPs that uses cutting-planes to approximate a convex feasible set.\n\n\n\n"
},

{
    "location": "library.html#Katana.KatanaSolver-Tuple{MathProgBase.SolverInterface.AbstractMathProgSolver}",
    "page": "Library",
    "title": "Katana.KatanaSolver",
    "category": "Method",
    "text": "KatanaSolver(lp_solver :: MathProgBase.AbstractMathProgSolver;\n             separator = KatanaFirstOrderSeparator(),\n             features = Vector{Symbol}(),\n             f_tol :: Float64 = 1e-6,\n             cut_coef_rng :: Float64 = 1e9,\n             log_level :: Int = 10,\n             iter_cap :: Int = 10000,\n             obj_eps :: Float64 = -1.0)\n\nConstruct a KatanaSolver with feasibility tolerance f_tol, a maximum coefficient range per cut cut_coef_rng and an iteration cap specifying the maximum number of rounds of LP solves + cut generation in iter_cap. Print out solver progress every log_level number of iterations, or suppress output with log_level=0. Use obj_eps as a stopping criterion: when the objective function evaluated by the LP has changed by less than obj_eps relative to the previous iteration.\n\ncut_coef_rng is used to round-off close-to-zero coefficients in generated cuts.\n\nThe separator is any implementing subtype of AbstractKatanaSeparator. It serves as the separation oracle used by the solver. The default is a first order separator that generates a single Newton cut per constraint.\n\nThe features vector is a list of optional features to enable in the solver. Currently supported are\n\n:VisData - Internal model logs actions to be exported and visualised\n\n\n\n"
},

{
    "location": "library.html#Katana.getKatanaCuts-Tuple{Katana.KatanaNonlinearModel}",
    "page": "Library",
    "title": "Katana.getKatanaCuts",
    "category": "Method",
    "text": "getKatanaCuts(m :: KatanaNonlinearModel)\n\nReturns a table of all linear cutting planes generated by Katana's separation oracle. If the model contains N variables (including an auxiliary variable) and m cutting planes were generated, the table will have size m x (N+2).\n\nFor every row: the first N columns represent the coefficients of model variables, the N+1th column a constant term, and the N+2th column an inequality direction (-1 for leq and 1 for geq). Thus, for coefficients a_1ldotsa_N, constant c and direction -1, re-construct the cutting plane as a_1x_1 + ldots + a_Nx_N leq c.\n\n\n\n"
},

{
    "location": "library.html#Katana.initialize!-Tuple{Katana.AbstractKatanaSeparator,Any,Any,Any,MathProgBase.SolverInterface.AbstractNLPEvaluator}",
    "page": "Library",
    "title": "Katana.initialize!",
    "category": "Method",
    "text": "initialize!(sep::AbstractKatanaSeparator, linear_model, num_var, num_constr, oracle::MathProgBase.AbstractNLPEvaluator)\n\nInitialise an instance of a subtype of AbstractKatanaSeparator with information about the KatanaNonlinearModel. This method is called by loadproblem! and MUST be overridden.\n\nlinear_model is the internal linear model of the KatanaNonlinearModel.\n\nnum_var is the number of solution variables, as passed by Katana. num_constr is the number of constraints in the problem, as passed by Katana. See solver implementation for details.\n\noracle can be queried for first- and second- derivative information and must be initialised in this method (see MathProgBase documentation on nonlinear models).\n\n\n\n"
},

{
    "location": "library.html#Katana.precompute!-Tuple{Katana.AbstractKatanaSeparator,Any}",
    "page": "Library",
    "title": "Katana.precompute!",
    "category": "Method",
    "text": "precompute!(sep::AbstractKatanaSeparator, xstar)\n\nImplement this method for a subtype of AbstractKatanaSeparator if your separator might only need to evaluate certain information once for all constraints using a solution from the internal LP model.\n\nxstar is the solution vector from the KatanaNonlinearModel's LP model\n\n\n\n"
},

{
    "location": "library.html#Katana.EpigraphNLPEvaluator",
    "page": "Library",
    "title": "Katana.EpigraphNLPEvaluator",
    "category": "Type",
    "text": "Wrapper for any MathProgBase.AbstractNLPEvaluator instance. Acts as an NLP evaluator for an NLP that has been transformed into epigraph form by treating the NL objective as if it were a constraint.\n\n\n\n"
},

{
    "location": "library.html#Katana.KatanaNonlinearModel",
    "page": "Library",
    "title": "Katana.KatanaNonlinearModel",
    "category": "Type",
    "text": "The MathProgBase non-linear model for Katana.\n\n\n\n"
},

{
    "location": "library.html#Katana.gencut-Tuple{Katana.AbstractKatanaSeparator,Any,Any,Any}",
    "page": "Library",
    "title": "Katana.gencut",
    "category": "Method",
    "text": "gencut!(sep::AbstractKatanaSeparator, xstar, i)\n\nGenerate a cut given an LP solution xstar for constraint i. bounds is a 2-tuple of (lb,ub). Since constraints are convex, one of the tuple bounds will be finite, and defines the level set of the constraint function. This method MUST be overridden for a subtype of AbstractKatanaSeparator. It is called by Katana as part of the solve routine.\n\nReturn a JuMP.AffExpr object.\n\n\n\n"
},

{
    "location": "library.html#Katana.isconstrsat-Tuple{Katana.AbstractKatanaSeparator,Any,Any,Any,Any}",
    "page": "Library",
    "title": "Katana.isconstrsat",
    "category": "Method",
    "text": "isconstrsat(sep::AbstractKatanaSeparator, i, lb, ub, f_tol)\n\nReturns true if the ith constraint is satisfied for the given bounds and tolerance. This method MUST be overriden for a subtype of AbstractKatanaSeparator as querying evaluated constraints is implementation-dependent.\n\n\n\n"
},

{
    "location": "library.html#Katana.numcuts-Tuple{Katana.KatanaNonlinearModel}",
    "page": "Library",
    "title": "Katana.numcuts",
    "category": "Method",
    "text": "numcuts(m::KatanaNonlinearModel)\n\nReturns the number of linear cuts added to the model. Includes linear constraints initially present.\n\n\n\n"
},

{
    "location": "library.html#Katana.numiters-Tuple{Katana.KatanaNonlinearModel}",
    "page": "Library",
    "title": "Katana.numiters",
    "category": "Method",
    "text": "numiters(m::KatanaNonlinearModel)\n\nReturns the number of iterations taken by the model.\n\n\n\n"
},

{
    "location": "library.html#Katana.jl-Library-1",
    "page": "Library",
    "title": "Katana.jl Library",
    "category": "section",
    "text": "Modules = [Katana]"
},

]}
