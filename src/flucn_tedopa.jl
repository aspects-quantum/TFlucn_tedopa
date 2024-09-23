module flucn_tedopa

export write_for_loop, write_to_file ##include("textfiles.jl")
export recur_coeff ##include("recur_coeffs.jl")
#export HBB_tf, HTOT_tf, exp_xHB_tdvp #include("v_ham_TDVP.jl") #  exp_xHB_tf, unit_gates_tf
#export exp_xHB_comp, tot_gate ###include("v_ham_TEBD.jl")

export exp_xHB_comp, tot_gate ###include("f_ham.jl")

using Reexport
@reexport using ITensors, ITensorMPS, PolyChaos, Plots, LinearAlgebra, QuadGK, LaTeXStrings, HDF5, Dates

include("textfiles.jl")
include("recur_coeffs.jl")
#include("v_ham_gates.jl")
#include("v_charfn.jl")
#include("v_thermal_MPS.jl")
#include("v_ham_TDVP.jl")
#include("v_ham_TEBD.jl")

include("f_ham.jl")
end