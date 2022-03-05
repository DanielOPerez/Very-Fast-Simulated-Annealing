include("vfsa_mod.jl")
include("cost_fun_mod.jl")
using .vfsa_mod
using .cost_fun_mod


npar=5
var1_op=Dict{Any,Any}("npar" => npar,
                      "xa" => -10.0,
                      "xb" => 10.0)

var2_op=Dict{Any,Any}("npar" => npar,
                      "xa" => -100.0,
                      "xb" => 100.0)

var_op=(var1_op,var2_op)
args_in=()

vfsa(cost_fun,
     var_op,
     itmax=30000,
     iop_prt=100,
     quench=1.0)


