using DrWatson
@quickactivate :flucn_tedopa
## https://juliadynamics.github.io/DrWatson.jl/stable/real_world/#Making-your-project-a-usable-module-1

#@disable_warn_order
ITensors.disable_warn_order()

#################            ########################################################################
################# Parameters ########################################################################
#################            ########################################################################



ρ = [0 0;0 1]
ρ = ρ/tr(ρ)

ITensors.op(::OpName"ρ", ::SiteType"S=1/2") = ρ


## time evolution parameters
cutoff = 1E-8
maxdim = 600
tau = 8*10^-3             ## time step duration
nt = 250
ttotal = nt*tau             ## TOTAL TIME evolution

N_chain = 120            ## Number of chain sites for single chain-transformed environment
tot_chain = 2*N_chain+1

boson_dim = 8

# Make an array of 'site' INDICES for the (spin+chain)
S_pos = N_chain+1

#s_total = [(n == S_pos_tilde) | (n == S_pos_real) ? Index(2, "S=1/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]
s_total = [(n == S_pos) ? Index(2, "S=1/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]


zer0 = zeros(boson_dim,boson_dim)
zer0[1,1] = 1
ITensors.op(::OpName"0", ::SiteType"Qudit", d::Int) = zer0


ω_C = 5                 ## bath cutoff
ω_0 = 0                 ## spin splitting
Ω = 1                   ## independent model if Ω = 0


u = 0.005 # counting field parameter


## 

p=plot()

t_list = collect(0:1:nt)*tau
T_list = [0.1]
α = 0.1

#cat(a, b) = reshape(append!(vec(a), vec(b)), size(a)[1:end-1]..., :)
support_cutoff = 6*10^2
supp = (0, support_cutoff)                 ## support of the weight function
Nquad = 10^7              ## Number of quadrature points
N_coeff = N_chain + 1

ab1 = Matrix{Float64}(undef,N_coeff,2)
ab2 = Matrix{Float64}(undef,N_coeff,2)


for T = T_list
    β = 1/T 
    n(k) = 1/(exp(β*k) - 1)
    w_fn1(k) = (2*α*k*exp(-k/ω_C))*(1 + n(k))           ## weight function for real space bath
    w_fn2(k) = (2*α*k*exp(-k/ω_C))*n(k)           ## weight function for tilde space bath
    η01 = quadgk(w_fn1, 0, support_cutoff)
    c_01 = sqrt(η01[1])
    η02 = quadgk(w_fn2, 0, support_cutoff)
    c_02 = sqrt(η02[1])
    ab1[1:90,1:2] = recur_coeff(w_fn1, supp, 90, Nquad)    ## recurrence coefficients for for real space bath
    ab2[1:90,1:2] = recur_coeff(w_fn2, supp, 90, Nquad)    ## recurrence coefficients for for tilde space bath
    if N_chain >= 90
        a1_100 = ab1[90,1]; b1_100 = ab1[90,2]; a2_100 = ab2[90,1]; b2_100 = ab2[90,2]
        [ab1[i,1] = a1_100  for i = 91:N_coeff]; [ab1[i,2] = b1_100  for i = 91:N_coeff]; 
        [ab2[i,1] = a2_100  for i = 91:N_coeff]; [ab2[i,2] = b2_100  for i = 91:N_coeff]; 
    end

    #chi_pu = charfn(ω_0, Ω, c_01, c_02, ab1, ab2, s_total, tau, nt, u, cutoff, maxdim)

    count_dates = exp_xHB(ab1, im*u, s_total);
    evol = unit_gates(ω_0, Ω, c_01, c_02, ab1, ab2, tau, s_total)

    ################################################################################################
    ρ = init_state(s_total)  ## initial state of (spin+chain);
    ################################################################################################

    chi_pu=Vector{ComplexF64}()

    U_ρ_Ud = ρ
    cd_U_ρ_Ud = apply(count_dates, U_ρ_Ud; cutoff, maxdim)

    chi = tr(cd_U_ρ_Ud)
        
    push!(chi_pu,chi[1])


    for t = 1:1:nt

        U_ρ_Ud = apply(evol, U_ρ_Ud; cutoff, maxdim, apply_dag=true)
        cd_U_ρ_Ud = apply(count_dates, U_ρ_Ud; cutoff, maxdim)

        chi = tr(cd_U_ρ_Ud)[1]

        @show imag(chi)/u

        
        push!(chi_pu,chi)

    end


    mean_Q = imag.(chi_pu)/u
    @show real(mean_Q)

    plot!((collect(0:1:length(mean_Q)-1)*tau)[begin:10:end], real(mean_Q)[begin:10:end], label= "T = $T, α = $α")

end


plot!(legend=:topright)
xlabel!("t")
title = string("<Q>, N_ch = ", N_chain,", b_dim = ", boson_dim,", u = ",u)
title!(title)
display("image/png", p)


a = "meanq_ub_1.png";
savefig(a)
    