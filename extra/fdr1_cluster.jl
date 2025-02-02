using DrWatson
@quickactivate :flucn_tedopa
## https://juliadynamics.github.io/DrWatson.jl/stable/real_world/#Making-your-project-a-usable-module-1

#@disable_warn_order
ITensors.disable_warn_order()

#################            ########################################################################
################# Parameters ########################################################################
#################            ########################################################################



ρ = [1 0; 0 0]
ρ = ρ / tr(ρ)

ITensors.op(::OpName"ρ", ::SiteType"S=1/2") = ρ


## time evolution parameters

cut = -11
cutoff = 10.0^cut
maxdim = 25
tau = 8*10^-3             ## time step duration
nt = 300
ttotal = nt * tau             ## TOTAL TIME evolution

N_chain = 150           ## Number of chain sites for single chain-transformed environment
tot_chain = 2*N_chain+1
S_pos = N_chain+1

println(S_pos)

n1_bsn_dim = 10;
b_dim = [n1_bsn_dim-round(Int64,(n1_bsn_dim-4.6)*(i-1)/(N_chain-1))  for i = 1:N_chain]   #       ## Dimension of chain sites
boson_dim = append!(reverse(b_dim),[0],b_dim)

s_total = [(n == S_pos) ? Index(2, "S=1/2") : Index(boson_dim[n], "Qudit") for n = 1:tot_chain]

ITensors.op(::OpName"0", ::SiteType"Qudit", d::Int) = 1.0I[1:d, 1] * 1.0I[1:d, 1]'
ITensors.op(::OpName"Idd", ::SiteType"Qudit", d::Int) = (1 / d) * Matrix(1.0I, d, d)
ITensors.op(::OpName"Idd", ::SiteType"S=1/2") = (1 / 2) * Matrix(1.0I, 2, 2)

state = [(n == S_pos) ? "ρ" : "0" for n = 1:tot_chain]
ρ = normalize(MPO(s_total, state))


ω_C = 5                 ## bath cutoff
ω_0 = 1                 ## spin splitting
Ω = 0                   ## independent model if Ω = 0
model = (ω_0 == 1) ? "independent" : "unbiased"

u = 0.005 # counting field parameter
t_list = collect(0:1:nt) * tau

# Reduce number of iterations for testing; increase later if needed.
T_list = [1, 4, 8]

α = 0.5
support_cutoff = 700
supp = (0, support_cutoff)                            ## support of the weight function
Nquad = 10^7                                          ## Reduced Number of quadrature points for speed
N_coeff = N_chain + 1

# Preallocate FDR to avoid dynamic resizing
FDR = Vector{Float64}(undef, length(T_list))

Threads.@threads for T_idx in 1:length(T_list)
    T = T_list[T_idx]
    ab1 = Matrix{Float64}(undef, N_coeff, 2)
    ab2 = Matrix{Float64}(undef, N_coeff, 2)

    β = 1 / T
    n(k) = 1 / (exp(β * k) - 1)
    w_fn1(k) = (2 * α * k * exp(-k / ω_C)) * (1 + n(k))            ## weight function for real space bath
    w_fn2(k) = (2 * α * k * exp(-k / ω_C)) * (n(k))                ## weight function for tilde space bath

    # Precompute coefficients to avoid recalculating them each time
    η001 = quadgk(w_fn1, 0, support_cutoff)
    c_01 = sqrt(Complex(η001[1]))
    η002 = quadgk(w_fn2, 0, support_cutoff)
    c_02 = sqrt(Complex(η002[1]))

    if N_chain >= 92
        ab1[1:92, 1:2] = recur_coeff(w_fn1, supp, 92, Nquad)    ## recurrence coefficients for for real space bath
        ab2[1:92, 1:2] = recur_coeff(w_fn2, supp, 92, Nquad)    ## recurrence coefficients for for tilde space bath
        a1_100, b1_100 = ab1[92, 1], ab1[92, 2]
        a2_100, b2_100 = ab2[92, 1], ab2[92, 2]
        fill!(ab1[93:end, 1], a1_100)
        fill!(ab1[93:end, 2], b1_100)
        fill!(ab2[93:end, 1], a2_100)
        fill!(ab2[93:end, 2], b2_100)
    else
        ab1 .= recur_coeff(w_fn1, supp, N_coeff, Nquad)    ## recurrence coefficients for real space bath
        ab2 .= recur_coeff(w_fn2, supp, N_coeff, Nquad)
    end

    count_gates = exp_xHB(ab1, ab2, im * u, s_total)
    evol = unit_gates(ω_0, Ω, c_01, c_02, ab1, ab2, tau, s_total)

    U_ρ_Ud = copy(ρ)  # avoid using global variable

    for t in 1:nt
        U_ρ_Ud = apply(evol, U_ρ_Ud; cutoff, maxdim, apply_dag=true)
        U_ρ_Ud = normalize(U_ρ_Ud)
        if t % 50 == 0  # Reduced frequency for calculating chi to save time
            println(t*tau)
            cd_U_ρ_Ud1 = apply(count_gates, U_ρ_Ud; cutoff, maxdim)
            chi_1 = tr(cd_U_ρ_Ud1)[1]
            mean_Q = real(imag(chi_1) / u)
            var_Q = real(-(log(conj(chi_1)) + log(chi_1)) / (u^2))
            @show mean_Q
            @show var_Q
        end
    end

    cd_U_ρ_Ud1 = apply(count_gates, U_ρ_Ud; cutoff, maxdim)
    chi_1 = tr(cd_U_ρ_Ud1)[1]
    mean_Q = real(imag(chi_1) / u)
    var_Q = real(-(log(conj(chi_1)) + log(chi_1)) / (u^2))

    FDR[T_idx] = var_Q / (T * mean_Q)

    @show FDR[T_idx]
end

@show FDR