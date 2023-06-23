using SOLNP
include("cost.jl")

# rosen
println("-----Rosen-----")
p0 = [0.0, 0.0]
res = solnp_solve(p=p0, cost=rosen)
println(res.p_out)

# cost
println("\n")
println("-----Cost-----")
p0 = Vector{Float64}([4.9, 0.1])
res = solnp_solve(p=p0, cost=cost)
println(res.p_out)