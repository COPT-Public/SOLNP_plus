include("SOLNP_CONST.jl")
using LinearAlgebra

_solnp_np::Int64 = 0
_solnp_nc::Int64 = 0
_solnp_nic::Int64 = 0
_solnp_nec::Int64 = 0
_solnp_cost_for_c = nothing
_solnp_grad_for_c = nothing
_solnp_hess_for_c = nothing

Base.@kwdef mutable struct SOLNPInput
    ibl::Vector{Float64}
    ibu::Vector{Float64}
    pbl::Vector{Float64}
    pbu::Vector{Float64}
    Ipc::Vector{Float64}
    Ipb::Vector{Float64}
    ib0::Vector{Float64}
    p::Vector{Float64}
    op::Vector{Float64}
    l::Vector{Float64}
    h::Vector{Float64}
    np::Int64
    nic::Int64
    nec::Int64
    nc::Int64
    cost
    grad
    hess
end

Base.@kwdef mutable struct SOLNPResult
    iter::Integer
    count_cost::Integer
    count_grad::Integer
    count_hess::Integer
    constraint::Float64
    restart_time::Integer
    obj::Float64
    status::Integer
    solve_time::Float64
    p_out::Vector{Float64}
    best_fea_p::Vector{Float64}
    ic::Vector{Float64}
    jh::Vector{Float64}
    ch::Vector{Float64}
    l_out::Vector{Float64}
    h_out::Vector{Float64}
end


function check_bound(; lb::Vector{Float64}, ub::Vector{Float64})::Bool
    reduce(&, lb .< ub)
end

function check_var_bound(; x::Vector{Float64}, lb::Vector{Float64}, ub::Vector{Float64})::Bool
    reduce(&, lb .< x) & reduce(&, x .< ub)
end

function cal_init_var(; lb::Vector{Float64}, ub::Vector{Float64})::Vector{Float64}
    x = (lb + ub) / 2
    x[isinf.(lb).&isinf.(ub)] .= 0.0
    x[isinf.(lb).&isfinite.(ub)] = ub[isinf.(lb).&isfinite.(ub)] .- 100
    x[isfinite.(lb).&isinf.(ub)] = lb[isfinite.(lb).&isinf.(ub)] .+ 100
    x
end

function checkAndFulfillProblem(m, ibl, ibu, ib0, pbl, pbu, p, cost, l, h)
    # init np, nic, nec, nc
    np = 0  # dim of variables
    nic = 0  # dim of inequality constraints
    nec = 0  # dim of equality constraints
    nc = 0  # dim of constraints

    # init Ipc
    Ipc = Vector{Float64}([0, 0])

    # init Ipb
    Ipb = Vector{Float64}([0, 0])

    # ibl
    if length(ibl) > 0
        nic = length(ibl)
    end

    # ibu
    if length(ibu) > 0
        nic = length(ibu)
    end

    # ib0
    if length(ib0) > 0
        nic = length(ib0)
    end

    # pbl
    if length(pbl) > 0
        np = length(pbl)
        Ipb[1] = 1.0
        Ipc[1] = 1.0
    end

    # pbu
    if length(pbu) > 0
        np = length(pbu)
        Ipb[1] = 1.0
        Ipc[2] = 1.0
    end

    # p
    if length(p) > 0
        np = length(p)
    end

    if Ipb[1] + nic >= 0.5
        Ipb[2] = 1.0
    end

    # check bounds
    if nic > 0
        if length(ibl) == 0
            ibl = Vector{Float64}(fill(-INFINITY, nic))
        end
        if length(ibu) == 0
            ibu = Vector{Float64}(fill(INFINITY, nic))
        end
        if !check_bound(lb=ibl, ub=ibu)
            error("ilb must be strictly less than ibu")
        elseif !check_var_bound(x=ib0, lb=ibl, ub=ibu)
            ib0 = cal_init_var(lb=ibl, ub=ibu)
        end
    end

    if np > 0
        if length(pbl) == 0
            pbl = Vector{Float64}(fill(-INFINITY, np))
        end
        if length(pbu) == 0
            pbu = Vector{Float64}(fill(INFINITY, np))
        end
        if !check_bound(lb=pbl, ub=pbu)
            error("pbl must be strictly less than pbu")
        elseif !check_var_bound(x=p, lb=pbl, ub=pbu)
            p = cal_init_var(lb=pbl, ub=pbu)
        end
    end

    # nc and nec
    if m >= 0
        nc = m
    else
        obs = cost(p)
        nc = length(obs) - 1
    end
    nec = nc - nic

    # l and h
    if length(l) != nc
        l = Vector{Float64}(fill(0.0, nc))
    end
    if length(h) != np^2
        h = Vector{Float64}(reshape(Matrix{Float64}(I, np, np), (:)))
    end

    # return
    return np, nic, nec, nc, Ipc, Ipb, ibl, ibu, ib0, pbl, pbu, p, l, h
end

function checkAndFulfillOption(option, np)::Vector{Float64}
    # default options
    default_opt = Vector([
        "rho" => 1.0,
        "pen_l1" => 1,
        "max_iter" => 50,
        "min_iter" => 10,
        "max_iter_rescue" => 50,
        "min_iter_rescue" => 10,
        "delta" => 1.0,
        "tol" => 1e-4,
        "tol_con" => 1e-3,
        "ls_time" => 10,
        "batchsize" => max(min(50, np ÷ 4), 1),
        "tol_restart" => 1.0,
        "re_time" => 5,
        "delta_end" => 1e-5,
        "maxfev" => 500 * np,
        "noise" => 1,
        "qpsolver" => 1,
        "scale" => 1,
        "bfgs" => 1,
        "rs" => 0,
        "grad" => 1,
        "k_i" => 3.0,
        "k_r" => 9,
        "c_r" => 10.0,
        "c_i" => 30.0,
        "ls_way" => 1,
        # 'rescue' => 1,
        "rescue" => 0,
        "drsom" => 0,
        "cen_diff":0,
        "gd_step":1e-1,
    ])

    # check input options
    for (i, kv) in enumerate(default_opt)
        if haskey(option, kv[1])
            default_opt[i] = kv[1] => option[kv[1]]
        end
    end

    return [kv[2] for kv in default_opt]
end

function wrapper_cost(ptr::Ptr{Float64}, n::Int64, result::Ptr{Float64}, n_res::Int64, par::Int32, _solnp_cost::Function)

    x = unsafe_wrap(Array, ptr, n)
    res_jl = unsafe_wrap(Array, result, n_res)
    c = _solnp_cost(x)

    if n_res > 1
        res_jl[:] = c
    else
        res_jl[1] = c
    end

    nothing
end

function wrapper_grad(ptr::Ptr{Float64}, n::Int64, result::Ptr{Float64}, n_res::Int64, _solnp_grad::Function)
    x = unsafe_wrap(Array, ptr, n)
    res_jl = unsafe_wrap(Array, result, n_res)
    g = _solnp_grad(x)
    res_jl[:] = reshape(g, (:))

    nothing
end

function wrapper_hess(ptr::Ptr{Float64}, n::Int64, result::Ptr{Float64}, n_res::Int64, _solnp_hess::Function)
    x = unsafe_wrap(Array, ptr, n)
    res_jl = unsafe_wrap(Array, result, n_res)
    h = _solnp_hess(x)
    res_jl[:] = reshape(h, (:))

    nothing
end

function checkAndWrapFunction(cost, grad, hess)
    # cost
    if cost !== nothing
        cost_for_c = @cfunction((ptr::Ptr{Float64}, result::Ptr{Float64}, par::Int32) -> wrapper_cost(ptr, _solnp_np, result, _solnp_nc + 1, par, _solnp_cost_for_c), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Cint))
    else
        error("cost must be a function")
    end

    # grad
    if grad !== nothing
        grad_for_c = @cfunction((ptr::Ptr{Float64}, result::Ptr{Float64}) -> wrapper_grad(ptr, _solnp_np, result, (_solnp_nc + 1) * _solnp_np, _solnp_cost_for_c), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}))
    else
        grad_for_c = Ptr{Cvoid}()
    end

    # hess
    if hess !== nothing
        hess_for_c = @cfunction((ptr::Ptr{Float64}, result::Ptr{Float64}) -> wrapper_hess(ptr, _solnp_np, result, (_solnp_nc + 1) * (_solnp_np^2), _solnp_cost_for_c), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}))
    else
        hess_for_c = Ptr{Cvoid}()
    end

    return cost_for_c, grad_for_c, hess_for_c
end

function set_global_vars(np_local, nc_local, nic_local, nec_local, cost_local, grad_local, hess_local)

    global _solnp_np = np_local
    global _solnp_nc = nc_local
    global _solnp_nic = nic_local
    global _solnp_nec = nec_local
    global _solnp_cost_for_c = cost_local
    global _solnp_grad_for_c = grad_local
    global _solnp_hess_for_c = hess_local

    nothing
end

function checkAndFulfillInputs(m, ibl, ibu, ib0, pbl, pbu, p, cost, grad, hess, option, l, h)::SOLNPInput
    # check and fulfill problems
    np, nic, nec, nc, Ipc, Ipb, ibl, ibu, ib0, pbl, pbu, p, l, h = checkAndFulfillProblem(m, ibl, ibu, ib0, pbl, pbu, p, cost, l, h)

    # set global variables
    set_global_vars(np, nc, nic, nec, cost, grad, hess)

    # check input options
    option = checkAndFulfillOption(option, np)

    # check input functions
    cost_for_c, grad_for_c, hess_for_c = checkAndWrapFunction(cost, grad, hess)


    # init solnp input
    solnp_input = SOLNPInput(ibl, ibu, pbl, pbu, Ipc, Ipb, ib0, p, option, l, h, np, nic, nec, nc, cost_for_c, grad_for_c, hess_for_c)
    solnp_input
end

function run_solnp_c(solnp_input::SOLNPInput)#::SOLNPResult
    # prepare for solve in c
    # output space
    scalars = zeros(Float64, 9)
    p_out = zeros(Float64, solnp_input.np)
    best_fea_p = zeros(Float64, solnp_input.np)
    ic = zeros(Float64, max(1, solnp_input.nic))
    jh = zeros(Float64, Int(solnp_input.op[OPT_INDEX["max_iter"]] + 1))
    ch = zeros(Float64, Int(solnp_input.op[OPT_INDEX["max_iter"]] + 1))
    l_out = zeros(Float64, max(1, solnp_input.nc))
    len_h = Bool(solnp_input.op[OPT_INDEX["bfgs"]]) ? (solnp_input.np + solnp_input.nic)^2 : 1
    h_out = zeros(Float64, len_h)

    # set cost
    # def_julia_callback_types = (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid})
    ccall((:def_python_callback, LIBSOLNP_PATH), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}), solnp_input.cost, solnp_input.grad, solnp_input.hess)

    # run solnp
    # solnp_c_types = (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble})

    ccall((:SOLNP_C, LIBSOLNP_PATH), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), solnp_input.ibl, solnp_input.ibu, solnp_input.pbl, solnp_input.pbu, solnp_input.Ipc, solnp_input.Ipb, solnp_input.ib0, solnp_input.p, solnp_input.op, solnp_input.l, solnp_input.h, solnp_input.np, solnp_input.nic, solnp_input.nec, scalars, p_out, best_fea_p, ic, jh, ch, l_out, h_out)

    # construct result
    solnp_result = SOLNPResult(Int(scalars[1]), Int(scalars[2]), Int(scalars[3]), Int(scalars[4]), scalars[5], Int(scalars[6]), scalars[7], Int(scalars[8]), scalars[9], p_out, best_fea_p, ic, jh, ch, l_out, h_out)
    solnp_result
end

function solnp_solve(; m::Int64=-1, ibl::Vector{Float64}=Vector{Float64}(), ibu::Vector{Float64}=Vector{Float64}(), ib0::Vector{Float64}=Vector{Float64}(), pbl::Vector{Float64}=Vector{Float64}(), pbu::Vector{Float64}=Vector{Float64}(), p::Vector{Float64}=Vector{Float64}(), cost, grad=nothing, hess=nothing, option::Dict{String,Float64}=Dict{String,Float64}(), l::Vector{Float64}=Vector{Float64}(), h::Matrix{Float64}=Matrix{Float64}(I, 0, 0))#::SOLNPResult
    # check and fulfill inputs
    solnp_input = checkAndFulfillInputs(m, ibl, ibu, ib0, pbl, pbu, p, cost, grad, hess, option, l, h)

    # solve problems
    solnp_result = run_solnp_c(solnp_input)

    # return results
    solnp_result
end
