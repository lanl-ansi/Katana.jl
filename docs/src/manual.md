# Katana User Manual

A more in-depth guide on customising, tweaking and extending Katana.

## Implementing Custom Separators

The Separator API is intended to abstract the logic of generating separating hyperplanes into two parts:
1. What stateful information is required? For example, the `KatanaFirstOrderSeparator` stores evaluated constraint
   values, evaluated constraint Jacobians, Jacobian sparsity structure, etc.
2. How is that information used to generate the cutting plane? In the default Newton-based algorithm, the hyperplane is
   constructed as a first-order expansion around a point.

Implementing a custom separator requires a distinction between the _information_ required and the _method_ used. It may not
be necessary to create a new subclass of `AbstractKatanaSeparator` if, for example, the separation oracle only intends to use
first-order information. Consider a method of generating cutting planes that involves starting at the LP optimal point ``x^{*}``
and performing gradient descent until within some ``\\epsilon`` of the constraint surface. Ideally, this would require another
function that can be used as the `algo` of a `KatanaFirstOrderSeparator`.

On the other hand, if the custom separator requires second-order constraint values, such as the Hessian, it would be necessary to implement a new subclass of `AbstractKatanaSeparator`. It is then up to that implementation if the actual algorithm called by `gencut` is itself modular.

## Visualising linear cuts

`KatanaNonlinearModel` provides a method to retrieve every linear cut generated during the solve process, as well as incremental LP solution vectors. Enable this feature by including `:VisData` in the list of features passed to `KatanaSolver`. Retrieve this information by calling
`getKatanaCuts` and `getKatanaSols` respectively. See the [library documentation](library.html#Katana.getKatanaCuts-Tuple{Katana.KatanaNonlinearModel}) for details.

An example approach to visualisation is provided [here](https://github.com/lanl-ansi/Katana.jl/blob/master/vis/vis.jl). Note that this script has external dependencies.

## Expressing models in JuMP for better performance with Katana

TODO.
