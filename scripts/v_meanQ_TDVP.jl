using DrWatson
@quickactivate :flucn_tedopa
## https://juliadynamics.github.io/DrWatson.jl/stable/real_world/#Making-your-project-a-usable-module-1

#@disable_warn_order
ITensors.disable_warn_order()

#################            ########################################################################
################# Parameters ########################################################################
#################            ########################################################################


ρ = [1 0;0 0]    
ρ = ρ/tr(ρ)             ## initial spin state

ITensors.op(::OpName"ρ", ::SiteType"S=1/2") = ρ


cutoff = 1E-9
#maxdim_ops = 5
maxdim = 100
tau = .01             ## time step duration
nt = 6
Nbeta = 10
ttotal = nt * tau           ## TOTAL TIME evolution

N_chain = 5            ## Number of chain sites for single chain-transformed environment
tot_chain = 2*N_chain+2


# Make an array of 'site' INDICES for the (spin+chain)
S_pos_t = N_chain+1
S_pos_r = N_chain+2

println(S_pos_r)

n1_bsn_dim = 7;
b_dim = [n1_bsn_dim-round(Int64,(n1_bsn_dim-2.6)*(i-1)/(N_chain-1)) for i = 1:N_chain]           ## Dimension of chain sites
boson_dim = append!(reverse(b_dim),[0 0],b_dim)

s_total = [(n == S_pos_t) | (n == S_pos_r) ? Index(2, "S=1/2") : Index(boson_dim[n], "Qudit") for n = 1:tot_chain]

#= function boson_zero(boson_dim)
    zer0 = zeros(boson_dim,boson_dim)
    zer0[1,1] = 1
    return zer0
end
ITensors.op(::OpName"0", ::SiteType"Qudit", d::Int) = boson_zero(d)
 =#
ITensors.op(::OpName"Idd", ::SiteType"Qudit", d::Int) = (1 / d) * Matrix(1.0I, d, d)
ITensors.op(::OpName"Idd", ::SiteType"S=1/2") = (1 / 2) * Matrix(1.0I, 2, 2)


ω_C = 5                 ## bath cutoff
ω_0 = 1                 ## spin splitting
Ω = 0                   ## independent model if Ω = 0

T = 5;                  ## temperature of bath
β = 1/T 
###

u = 0.01 # counting field parameter


## 

p=plot()

t_list = collect(0:1:nt)*tau
#α_list = [0.1,1.5]
α = .5

#cat(a, b) = reshape(append!(vec(a), vec(b)), size(a)[1:end-1]..., :)
support_cutoff = 7*10^2
supp = (0, support_cutoff)                            ## support of the weight function
Nquad = 10^7                                          ## Number of quadrature points
N_coeff = N_chain + 1

ab1 = Matrix{Float64}(undef,N_coeff,2)
ab2 = Matrix{Float64}(undef,N_coeff,2)

n(k) = 1/(exp(β*k) - 1)

#for α = α_list
    w_fn1(k) = (2*α*k*exp(-k/ω_C))*(1 + n(k))           ## weight function for real space bath
    w_fn2(k) = (2*α*k*exp(-k/ω_C))*n(k)           ## weight function for tilde space bath
    η01 = quadgk(w_fn1, 0, support_cutoff)
    c_01 = sqrt(η01[1])
    η02 = quadgk(w_fn2, 0, support_cutoff)
    c_02 = sqrt(η02[1])
    
    if N_chain >= 90
        ab1[1:90,1:2] = recur_coeff(w_fn1, supp, 90, Nquad)    ## recurrence coefficients for for real space bath
        ab2[1:90,1:2] = recur_coeff(w_fn2, supp, 90, Nquad)    ## recurrence coefficients for for tilde space bath
        a1_100 = ab1[90,1]; b1_100 = ab1[90,2]; a2_100 = ab2[90,1]; b2_100 = ab2[90,2]
        [ab1[i,1] = a1_100  for i = 91:N_coeff]; [ab1[i,2] = b1_100  for i = 91:N_coeff]; 
        [ab2[i,1] = a2_100  for i = 91:N_coeff]; [ab2[i,2] = b2_100  for i = 91:N_coeff]; 
    else
        ab1 = recur_coeff(w_fn1, supp, N_coeff, Nquad)    ## recurrence coefficients for for real space bath
        ab2 = recur_coeff(w_fn2, supp, N_coeff, Nquad)    ## recurrence coefficients for for tilde space bath
    end

    bath_ham = HB_tf(ab1, s_total)
    tot_ham = HTOT_tf(ω_0,Ω,c_01,c_02,ab1, ab2, s_total)

    ################################################################################################
    ρ_vec = thermal_MPS(ab1, s_total, β, Nbeta, cutoff, maxdim)  ## initial state of (spin+chain);
    ################################################################################################

    char_fn=Vector{ComplexF64}()

    cg_ρ = tdvp(bath_ham, -u*1.0im, ρ_vec; nsteps=5, updater_backend="exponentiate", cutoff=1e-9, maxdim, normalize=true)
    U_cg_ρ = cg_ρ
    U_ρ = ρ_vec
    cg_U_ρ = cg_ρ

    chi = inner(cg_U_ρ, U_cg_ρ)
    println(chi[1])
    push!(char_fn,chi[1])

    for t in 1:nt
        U_cg_ρ = tdvp(tot_ham, -1.0im*tau, U_cg_ρ; nsteps=10, updater_backend="exponentiate", cutoff=1e-9, maxdim, normalize=true)
        U_ρ = tdvp(tot_ham, -1.0im*tau, U_ρ; nsteps=10, updater_backend="exponentiate", cutoff=1e-9, maxdim, normalize=true)
        cg_U_ρ = tdvp(bath_ham, -1.0im*u, U_ρ; nsteps=5, updater_backend="exponentiate", cutoff=1e-9, maxdim, normalize=true)

        chi = inner(cg_U_ρ, U_cg_ρ)
        println(chi[1])
        push!(char_fn,chi[1])
    end
chi_pu = char_fn
    #chi_pu = charfn_TDVP(ω_0, Ω, c_01, c_02, ab1, ab2, s_total, β, Nbeta, tau, nt, u, cutoff, maxdim)

    mean_Q = imag.(chi_pu)/u
    @show real(mean_Q)

    #plot!(t_list,real(mean_Q),label= "α = $α")

    plot!((collect(0:1:length(mean_Q)-1)*tau),real(mean_Q),label= "α = $α")
#end


plot!(legend=:topright)
xlabel!("t")
title = string("<Q>, N_ch = ", N_chain,", b_dim = ", boson_dim,", u = ",u, ", k_max = ", support_cutoff)
title!(title)
display("image/png", p)

#= 
a = "meanq_tdvp_1.png";
savefig(a)
     =#
#= 
using HDF5
f = h5open("myfile.h5","r")
T = read(f,"thermal_state",MPS)
close(f) =#