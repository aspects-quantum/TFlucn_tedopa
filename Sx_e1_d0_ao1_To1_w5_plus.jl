using DrWatson
@quickactivate :flucn_tedopa

ITensors.disable_warn_order()

##########################################################################
function dw_unit_gates(ω_0, Ω, c_01, c_02, ab1, ab2, tau, s_total)
        tot_chain = length(s_total)
        S_pos = Int((tot_chain + 1) / 2)

        ω_n_REAL = ab1[1:S_pos-1, 1]
        ω_n_TILD = ab2[1:S_pos-1, 1]
        ω_n_total = append!(reverse(ω_n_TILD), [0], ω_n_REAL)
        t_n_REAL = sqrt.(ab1[2:S_pos, 2])
        t_n_TILD = sqrt.(ab2[2:S_pos, 2])
        t_n_total = append!(reverse(t_n_TILD), [0], t_n_REAL)

        gates = ITensor[]

        for j in 2:tot_chain-1

                if j < S_pos
                        s1 = s_total[j]
                        s2 = s_total[j-1]

                        ω_n = ω_n_total[j]
                        t_n = t_n_total[j]
                        hj = (-ω_n) * op("N", s1) * op("Id", s2) +
                                 (-t_n) * op("Adag", s1) * op("A", s2) +
                                 (-t_n) * op("A", s1) * op("Adag", s2)
                        Gj = exp(-im * (tau / 2) * hj)
                        push!(gates, Gj)

                elseif j == S_pos
                        hj = ω_0 * op("Sz", s_total[j]) * op("Id", s_total[j+1]) +
                                 Ω * op("Sx", s_total[j]) * op("Id", s_total[j+1]) +
                                 c_01 * op("Sx", s_total[j]) * op("A", s_total[j+1]) +
                                 c_02 * op("Sx", s_total[j]) * op("Adag", s_total[j+1])
                        Gj = exp(-im * (tau / 2) * hj)
                        push!(gates, Gj)
                        hj = c_02 * op("Sx", s_total[j]) * op("A", s_total[j-1]) +
                                 c_02 * op("Sx", s_total[j]) * op("Adag", s_total[j-1])
                        Gj = exp(-im * (tau / 2) * hj)
                        push!(gates, Gj)

                else
                        s1 = s_total[j]
                        s2 = s_total[j+1]

                        ω_n = ω_n_total[j]
                        t_n = t_n_total[j]
                        hj = ω_n * op("N", s1) * op("Id", s2) +
                                 t_n * op("Adag", s1) * op("A", s2) +
                                 t_n * op("A", s1) * op("Adag", s2)
                        Gj = exp(-im * (tau / 2) * hj)
                        push!(gates, Gj)
                end
        end

        return append!(gates, reverse(gates))

end


# Method definitions must be at the top level, not inside functions
ITensors.op(::OpName"ρ", ::SiteType"S=1/2") = [1. 1.; 1. 1.] ./2 # Adjusted normalization of the spin state
ITensors.op(::OpName"0", ::SiteType"Qudit", d::Int) = 1.0I[1:d, 1] * 1.0I[1:d, 1]'
ITensors.op(::OpName"Idd", ::SiteType"Qudit", d::Int) = (1 / d) * Matrix(1.0I, d, d)
ITensors.op(::OpName"Idd", ::SiteType"S=1/2") = (1 / 2) * Matrix(1.0I, 2, 2)

let
        ################# Parameters ########################################################################
        # Define filename for output
        file_name_txt_obs = string(split(split(@__FILE__, ".")[end-1], string('\\'))[end], "_obs.txt")
        @show file_name_txt_obs

        # Define parameters for simulation
        cut = -11  # Cutoff for singular values
        cutoff = 10.0^cut
        maxdim = 80
        tau = 0.001  # Time step duration
        jump = 20
        nt = 10000  # Number of time steps
        ttotal = nt * tau  # Total time evolution

        N_chain = 180  # Number of chain sites for a single chain-transformed environment
        tot_chain = 2 * N_chain + 1
        S_pos = N_chain + 1

        println(S_pos)

        n1_bsn_dim = 14  # Dimension of chain sites
        b_dim_real = [n1_bsn_dim - round(Int64, (n1_bsn_dim - 3.6) * (i - 1) / (N_chain - 1)) for i in 1:N_chain]  # Dimension of chain sites
        b_dim_tilde = [Int(n1_bsn_dim - 4 - round(Int64, (n1_bsn_dim - 4 - 1.6) * (i - 1) / (N_chain - 1))) for i in 1:N_chain]
        boson_dim = append!(reverse(b_dim_tilde), [0], b_dim_real)

        s_total = [(n == S_pos) ? Index(2, "S=1/2") : Index(boson_dim[n], "Qudit") for n in 1:tot_chain]

        state = [(n == S_pos) ? "ρ" : "0" for n in 1:tot_chain]
        ρ = normalize(MPO(s_total, state))

        # Bath parameters
        ω_C = 5  # Bath cutoff
        ω_0 = 1  # Spin splitting
        Ω = 0
        model = (ω_0 == 1) ? "local" : "tunnel"
        model = (Ω == 1) ? "tunnel" : "local"

        T = .1  # Temperature of bath
        β = 1 / T

        t_list = collect(0:1:nt) * tau
        α = 0.1

        support_cutoff = 700
        supp = (0, support_cutoff)  # Support of the weight function
        Nquad = 10^7  # Reduced number of quadrature points for speed
        N_coeff = N_chain + 1
        N_rec = 92
        ab1 = Matrix{Float64}(undef, N_coeff, 2)
        ab2 = Matrix{Float64}(undef, N_coeff, 2)

        # Define functions for the weight functions
        n(ω) = 1 / (exp(β * ω) - 1)
        w_fn1(k) = (2 * α * k * exp(-k / ω_C)) * (1 + n(k))
        w_fn2(k) = (2 * α * k * exp(-k / ω_C)) * n(k)

        # Calculate recurrence coefficients
        η01 = quadgk(w_fn1, 0, support_cutoff)
        c_01 = sqrt(Complex(η01[1]))
        η02 = quadgk(w_fn2, 0, support_cutoff)
        c_02 = sqrt(Complex(η02[1]))

        if N_chain >= N_rec
                ab1[1:N_rec, 1:2] = recur_coeff(w_fn1, supp, N_rec, Nquad)
                ab2[1:N_rec, 1:2] = recur_coeff(w_fn2, supp, N_rec, Nquad)
                a1_100, b1_100 = ab1[N_rec, 1], ab1[N_rec, 2]
                a2_100, b2_100 = ab2[N_rec, 1], ab2[N_rec, 2]
                ab1[N_rec+1:N_coeff, 1:2] .= repeat([a1_100 b1_100], N_coeff - N_rec, 1)
                ab2[N_rec+1:N_coeff, 1:2] .= repeat([a2_100 b2_100], N_coeff - N_rec, 1)
        else
                ab1 .= recur_coeff(w_fn1, supp, N_coeff, Nquad)
                ab2 .= recur_coeff(w_fn2, supp, N_coeff, Nquad)
        end

        OBS = MPO(OpSum() + (1, "Sx", S_pos), s_total)
        evol = dw_unit_gates(ω_0, Ω, c_01, c_02, ab1, ab2, tau, s_total)

        # Initialize characteristic function vector
        t_plot = Float64[]

        obs_list = Float64[]

        # Time evolution of density matrix
        U_ρ_Ud = ρ
        O_U_ρ_Ud = apply(OBS, U_ρ_Ud; cutoff = 0.1*cutoff)
        @show obs = real(tr(O_U_ρ_Ud))
        push!(obs_list, obs)
        write_for_loop(file_name_txt_obs, string(1), "$(model) boson: T = $T, alpha = $α, N_chain = $N_chain, maxdim = $maxdim, cutoff = $cut, tau = $tau, jump = $jump, omega_C = $ω_C, boson_dim = $n1_bsn_dim")
        write_for_loop(file_name_txt_obs, string(2), string(obs))

        for t in 1:nt
                U_ρ_Ud = normalize(apply(evol, U_ρ_Ud; cutoff, maxdim, apply_dag = true))
                #U_ρ_Ud = add(0.5 * swapprime(dag(U_ρ_Ud), 0 => 1), 0.5 * U_ρ_Ud; maxdim = 2 * maxdim)
                #U_ρ_Ud = normalize(ITensors.truncate(ITensors.truncate(U_ρ_Ud; maxdim = 5, site_range = (1:(Int(3*floor(N_chain / 4))))); maxdim = 20, site_range = (tot_chain:tot_chain-(Int(floor(3*(N_chain / 4)))))))
                orthogonalize!(U_ρ_Ud, S_pos)
                if t % jump == 1
                        U_ρ_Ud /= tr(U_ρ_Ud)
                        O_U_ρ_Ud = apply(OBS, U_ρ_Ud; cutoff = cutoff)

                        @show obs = real(tr(O_U_ρ_Ud))
                        write_for_loop(file_name_txt_obs, string(t + 1), string(obs))
                        push!(obs_list, obs)
                end

                @show t * tau
                @show maxlinkdim(U_ρ_Ud)
        end

        write_to_file(file_name_txt_obs, "$(model) boson: T = $T, alpha = $α, N_chain = $N_chain, maxdim = $maxdim, cutoff = $cut, tau = $tau, omega_C = $ω_C, boson_dim = $n1_bsn_dim", string(obs_list))
end