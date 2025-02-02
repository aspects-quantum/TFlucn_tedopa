##########################################################################################################
############# calculate MEAN HEAT EXCHANGED with the bosonic bath at any time t  #########################
##########################################################################################################



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


## time evolution parameters
cutoff = 1E-8
maxdim = 1000
tau = 8*10^-3             ## time step duration
nt = 200
ttotal = nt*tau           ## TOTAL TIME evolution

N_chain = 120             ## Number of chain sites for single chain-transformed environment
tot_chain = 2*N_chain+1 

boson_dim = 8            ## boson excitation cutoff

# Make an array of 'site' INDICES for the (spin+chain)
S_pos = N_chain+1
s_total = [(n == S_pos) ? Index(2, "S=1/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]

    
zer0 = zeros(boson_dim,boson_dim)
zer0[1,1] = 1
ITensors.op(::OpName"0", ::SiteType"Qudit", d::Int) = zer0


ω_C = 5                   ## bath cutoff
ω_0 = 1                   ## spin splitting
Ω = 0                     ## independent model if Ω = 0, unbiased model if Ω = 1

T = .1;                   ## temperature of bath
β = 1/T 
###

u = 0.01                                              ## counting field parameter

t_list = collect(0:1:nt)*tau
#α_list = [0.1,1.5]
α = .2

support_cutoff = 10^4
supp = (0, support_cutoff)                            ## support of the weight function
Nquad = 10^7                                          ## Number of quadrature points
N_coeff = N_chain + 1

ab1 = Matrix{Float64}(undef,N_coeff,2)
ab2 = Matrix{Float64}(undef,N_coeff,2)

n(k) = 1/(exp(β*k) - 1)

p = plot()

#for α = α_list
    w_fn1(k) = (2*α*k*exp(-k/ω_C))*(1 + n(k))         ## weight function for real space bath
    w_fn2(k) = (2*α*k*exp(-k/ω_C))*n(k)               ## weight function for tilde space bath
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

    char_fn=Vector{ComplexF64}()

    U_ρ_Ud = ρ
    cd_U_ρ_Ud = apply(count_dates, U_ρ_Ud; cutoff, maxdim)

    chi = tr(cd_U_ρ_Ud)
        
    push!(char_fn,chi[1])


    for t in 1:nt
    
        global U_ρ_Ud = apply(evol, U_ρ_Ud; cutoff, maxdim, apply_dag=true)
        global cd_U_ρ_Ud = apply(count_dates, U_ρ_Ud; cutoff, maxdim)

        @show chi = tr(cd_U_ρ_Ud)
        
        push!(char_fn,chi[1])

    end
    chipu = char_fn
    mean_Q = imag.(chi_pu)/u
    @show real(mean_Q)

    #plot!((collect(0:1:length(mean_Q)-1)*tau)[begin:5:end],real(mean_Q)[begin:5:end],label= "α = $α")
    plot!((collect(0:1:length(mean_Q)-1)*tau),real(mean_Q),label= "α = $α")
#end


plot!(legend=:topright)
xlabel!("t")
title = string("<Q>, N_ch = ", N_chain,", b_dim = ", boson_dim,", u = ",u, ", k_max = ", support_cutoff)
title!(title)
display("image/png", p)

#= 
a = "meanq_.png";
savefig(a)
    
 =#