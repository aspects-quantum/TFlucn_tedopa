using DrWatson
@quickactivate :flucn_tedopa

ITensors.disable_warn_order()

# Method definitions must be at the top level, not inside functions
ITensors.op(::OpName"ρ", ::SiteType"S=1/2") = [1 0; 0 0]  # Adjusted normalization of the spin state
ITensors.op(::OpName"0", ::SiteType"Qudit", d::Int) = 1.0I[1:d, 1] * 1.0I[1:d, 1]'
ITensors.op(::OpName"Idd", ::SiteType"Qudit", d::Int) = (1 / d) * Matrix(1.0I, d, d)
ITensors.op(::OpName"Idd", ::SiteType"S=1/2") = (1 / 2) * Matrix(1.0I, 2, 2)

################# Parameters ########################################################################
# Define filename for output
file_name_txt = string(split(split(@__FILE__, ".")[end-1], string('\\'))[end], ".txt")
@show file_name_txt

# Define parameters for simulation
cut = -11
cutoff = 10.0^cut
maxdim = 40

N_chain = 120  # Number of chain sites for a single chain-transformed environment
tot_chain = 2 * N_chain + 1
S_pos = N_chain + 1

println(S_pos)

n1_bsn_dim = 11
b_dim_real = [n1_bsn_dim - round(Int64, (n1_bsn_dim - 3.6) * (i - 1) / (N_chain - 1)) for i in 1:N_chain]  # Dimension of chain sites
b_dim_tilde = [Int(n1_bsn_dim - 2 - round(Int64, (n1_bsn_dim - 2 - 1.6) * (i - 1) / (N_chain - 1))) for i in 1:N_chain]
boson_dim = append!(reverse(b_dim_tilde), [0], b_dim_real)

s_total = [(n == S_pos) ? Index(2, "S=1/2") : Index(boson_dim[n], "Qudit") for n in 1:tot_chain]

state = [(n == S_pos) ? "ρ" : "0" for n in 1:tot_chain]
ρ = normalize(MPO(s_total, state))

# Bath parameters
ω_C = 5  # Bath cutoff
ω_0 = 0  # Spin splitting
Ω = 1  # Independent model if Ω = 0
model = (ω_0 == 1) ? "independent" : "unbiased"

# Reduce number of iterations for testing; increase later if needed.
T_list = [3, 5, 8, 10]

α = 1
support_cutoff = 700
supp = (0, support_cutoff)                            ## support of the weight function
Nquad = 10^7                                          ## Reduced Number of quadrature points for speed
N_coeff = N_chain + 1

# Preallocate FDR to avoid dynamic resizing
FDR = Vector{Float64}(undef, length(T_list))
mQ = Vector{Float64}(undef, length(T_list))
vQ = Vector{Float64}(undef, length(T_list))

Threads.@threads for T_idx in 1:length(T_list)
        T = T_list[T_idx]
        ab1 = Matrix{Float64}(undef, N_coeff, 2)
        ab2 = Matrix{Float64}(undef, N_coeff, 2)

        tau = 0.002*10/T  # Time step duration
        u = 2*tau

        β = 1 / T
        n(k) = 1 / (exp(β * k) - 1)
        w_fn1(k) = (2 * α * k * exp(-k / ω_C)) * (1 + n(k))            ## weight function for real space bath
        w_fn2(k) = (2 * α * k * exp(-k / ω_C)) * (n(k))                ## weight function for tilde space bath

        # Calculate recurrence coefficients
        η01 = quadgk(w_fn1, 0, support_cutoff)
        c_01 = sqrt(Complex(η01[1]))
        η02 = quadgk(w_fn2, 0, support_cutoff)
        c_02 = sqrt(Complex(η02[1]))

        if N_chain >= 92
                ab1[1:92, 1:2] = recur_coeff(w_fn1, supp, 92, Nquad)
                ab2[1:92, 1:2] = recur_coeff(w_fn2, supp, 92, Nquad)
                a1_100, b1_100 = ab1[92, 1], ab1[92, 2]
                a2_100, b2_100 = ab2[92, 1], ab2[92, 2]
                ab1[93:N_coeff, 1:2] .= repeat([a1_100 b1_100], N_coeff - 92, 1)
                ab2[93:N_coeff, 1:2] .= repeat([a2_100 b2_100], N_coeff - 92, 1)
        else
                ab1 .= recur_coeff(w_fn1, supp, N_coeff, Nquad)
                ab2 .= recur_coeff(w_fn2, supp, N_coeff, Nquad)
        end

        # Define counting gates and evolution gates
        count_gates = exp_xHB(ab1, ab2, im * u, s_total)
        evol = unit_gates(ω_0, Ω, c_01, c_02, ab1, ab2, tau, s_total)

        # Initialize characteristic function vector
        char_fn = Vector{ComplexF64}()

        gamma = (pi/4) * (2*α*Ω * exp(-Ω / ω_C)) * coth(β*Ω/2)
        write_for_loop(file_name_txt, string(1), "t_inf=$(5/gamma)")
        @show 5/gamma

        # Time evolution of density matrix
        U_ρ_Ud = copy(ρ)
        write_for_loop(file_name_txt, string(1), "$(model) boson: T = $T, alpha = $α, N_chain = $N_chain, maxdim = $maxdim, cutoff = $cut, tau = $tau, u = $u, boson_dim = $n1_bsn_dim")

        t = 0
        i = 1
        while t < 1
                U_ρ_Ud = apply(evol, U_ρ_Ud; cutoff, maxdim, apply_dag = true)
                if i%10 == 0
                        @show t
                        U_ρ_Ud = add(0.5 * swapprime(dag(U_ρ_Ud), 0 => 1), 0.5 * U_ρ_Ud; maxdim = 300)
                        U_ρ_Ud = normalize(ITensors.truncate(ITensors.truncate(U_ρ_Ud; maxdim = 10, site_range = (1:(Int(floor(3*N_chain / 4))))); maxdim = 15, site_range = (tot_chain:tot_chain-(Int(floor(3*(N_chain / 4)))))))
                end

                orthogonalize!(U_ρ_Ud, findfirst(==(maxlinkdim(U_ρ_Ud)), linkdims(U_ρ_Ud)))
                t += tau
                i += 1
        end
        @show t
        cd_U_ρ_Ud = apply(count_gates, U_ρ_Ud; cutoff = cutoff * 0.1)
        chi = tr(cd_U_ρ_Ud)[1]
        @show mQ[T_idx] = real(imag(chi) / u)
        @show vQ[T_idx] = real(-(log(conj(chi)) + log(chi)) / (u^2))
        @show FDR[T_idx] = vQ[T_idx] / (T * mQ[T_idx])
        write_for_loop(file_name_txt, string(2), "T = $T: mQ = $(mQ[T_idx]), var_Q = $(vQ[T_idx]), FDR = $(FDR[T_idx])")
end
@show FDR
write_to_file(file_name_txt, "$(model) boson: T = $(T_list), alpha = $α, N_chain = $N_chain, maxdim = $maxdim, cutoff = $cut, tau = $tau, u = $u, boson_dim = $n1_bsn_dim", string(real(FDR)))
write_to_file(file_name_txt, "$(model) boson: T = $(T_list), alpha = $α, N_chain = $N_chain, maxdim = $maxdim, cutoff = $cut, tau = $tau, u = $u, boson_dim = $n1_bsn_dim", string(real(vQ)))
write_to_file(file_name_txt, "$(model) boson: T = $(T_list), alpha = $α, N_chain = $N_chain, maxdim = $maxdim, cutoff = $cut, tau = $tau, u = $u, boson_dim = $n1_bsn_dim", string(real(mQ)))