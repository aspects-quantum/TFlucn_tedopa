using DrWatson
@quickactivate :flucn_tedopa
using ITensors

ITensors.disable_warn_order()

#################### Parameters ####################
# Initial density matrix
ρ_init = [1 0; 0 0]
ρ_init /= tr(ρ_init)

# Define the density matrix operator for ITensors
ITensors.op(::OpName"ρ", ::SiteType"S=1/2") = ρ_init

# Time evolution parameters
cut_exponent = -11
cutoff = 10.0^cut_exponent
maxdim = 100
tau = .01             # Time step duration
nt = 300               # Number of time steps
ttotal = nt * tau      # Total time evolution

# Chain and bosonic site parameters
N_chain = 90                          # Number of chain sites
tot_chain = 2 * N_chain + 1
S_pos = N_chain + 1                   # Position of spin in the chain
println("Spin position in chain: ", S_pos)

n1_bsn_dim = 10
b_dim = [n1_bsn_dim - round(Int, (n1_bsn_dim - 4.6) * (i - 1) / (N_chain - 1)) for i in 1:N_chain]
boson_dim = append!(reverse(b_dim), [0], b_dim)

# Site indices for tensor network
s_total = [(n == S_pos) ? Index(2, "S=1/2") : Index(boson_dim[n], "Qudit") for n in 1:tot_chain]

# Bosonic operators
ITensors.op(::OpName"0", ::SiteType"Qudit", d::Int) = 1.0I[1:d, 1] * 1.0I[1:d, 1]'
ITensors.op(::OpName"Idd", ::SiteType"Qudit", d::Int) = (1 / d) * Matrix(1.0I, d, d)
ITensors.op(::OpName"Idd", ::SiteType"S=1/2") = (1 / 2) * Matrix(1.0I, 2, 2)

# Model parameters
ω_C = 5.0        # Bath cutoff frequency
ω_0 = 0.0        # Spin splitting
Ω = 1.0          # Model parameter
model = (ω_0 == 1) ? "independent" : "unbiased"
u = 0.005       # Counting field parameter (for variance)

# Time and thermal parameters
t_list = collect(0:1:nt) * tau
α = 1.0                          # Coupling strength
support_cutoff = 700.0           # Support for weight function
Nquad = 10^7                     # Number of quadrature points
N_coeff = N_chain + 1            # Number of recurrence coefficients
T_list = [1.0, 4.0]              # Temperature list

# Initial state and normalization
state = [(n == S_pos) ? "ρ" : "0" for n in 1:tot_chain]
ρ = normalize(MPO(s_total, state))

# Function to compute thermal occupation number
n(k, β) = 1.0 / (exp(β * k) - 1.0)

# Initialize data collection
FDR = Vector{Float64}()

#################### Main Calculation Loop ####################
for T in T_list
    β = 1.0 / T

    # Define weight functions for bath
    w_fn1(k) = (2 * α * k * exp(-k / ω_C)) * (1.0 + n(k, β))   # Real-space bath weight function
    w_fn2(k) = (2 * α * k * exp(-k / ω_C)) * n(k, β)           # Tilde-space bath weight function

    # Calculate quadrature and recurrence coefficients
    η01 = quadgk(w_fn1, 0.0, support_cutoff)[1]
    c_01 = sqrt(η01)
    η02 = quadgk(w_fn2, 0.0, support_cutoff)[1]
    c_02 = sqrt(η02)

    # Initialize recurrence coefficients
    ab1 = recur_coeff(w_fn1, (0.0, support_cutoff), N_coeff, Nquad)
    ab2 = recur_coeff(w_fn2, (0.0, support_cutoff), N_coeff, Nquad)

    # Set up counting gates and evolution operator
    count_gates = exp_xHB(ab1, ab2, im * u, s_total)
    evol = unit_gates(ω_0, Ω, c_01, c_02, ab1, ab2, tau, s_total)

    # Calculate characteristic function at initial time
    U_ρ_Ud = ρ
    cd_U_ρ_Ud = apply(count_gates, ρ; cutoff, maxdim)
    chi_1 = tr(cd_U_ρ_Ud)
    chi_2 = conj(chi_1)
    variance = real(-(log(chi_2[1]) + log(chi_1[1])) / (u^2))
    println("Variance at T = $T, initial time: ", variance)

    # Store initial variance
    push!(FDR, variance)

    # Time evolution loop
    for t in 1:nt
        # Apply time evolution
        U_ρ_Ud = apply(evol, U_ρ_Ud; cutoff, maxdim, apply_dag=true)
        U_ρ_Ud = normalize(U_ρ_Ud)

        # Update characteristic function and variance at each time step
        cd_U_ρ_Ud = apply(count_gates, U_ρ_Ud; cutoff, maxdim)
        chi_1 = tr(cd_U_ρ_Ud)
        chi_2 = conj(chi_1)
        variance = real(-(log(chi_2[1]) + log(chi_1[1])) / (u^2))
        println("Variance at T = $T, time = $(t * tau): ", variance)
    end
end
