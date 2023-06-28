LIBSOLNP_PATH = "libsolnp"

# INFINITY = 1e20
INFINITY = Inf64
TOL = 1e-8

OPT_INDEX = Dict{String,Int}([
    "rho" => 1,
    "pen_l1" => 2,
    "max_iter" => 3,
    "min_iter" => 4,
    "max_iter_rescue" => 5,
    "min_iter_rescue" => 6,
    "delta" => 7,
    "tol" => 8,
    "tol_con" => 9,
    "ls_time" => 10,
    "batchsize" => 11,
    "tol_restart" => 12,
    "re_time" => 13,
    "delta_end" => 14,
    "maxfev" => 15,
    "noise" => 16,
    "qpsolver" => 17,
    "scale" => 18,
    "bfgs" => 19,
    "rs" => 20,
    "grad" => 21,
    "k_i" => 22,
    "k_r" => 23,
    "c_r" => 24,
    "c_i" => 25,
    "ls_way" => 26,
    # 'rescue' => 1,
    "rescue" => 27,
    "drsom" => 28,
    "cen_diff": 29,
    "gd_step": 30,
])