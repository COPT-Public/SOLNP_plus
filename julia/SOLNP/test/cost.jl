"""
    cost function
        rosenbrock
"""
function rosen(x::Vector)
    return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

"""
    cost function
```math
    f(x) = (x_1 - 5)^2 + x_2^2 - 25
    g(x) = -x_1^2 + x_1
```
"""
function cost(x::Array{Float64})
    # x is a vector of length 2
    f1 = (x[1] - 5)^2 + x[2]^2 - 25
    f2 = -x[1]^2 + x[1]
    [f1, f2]
end