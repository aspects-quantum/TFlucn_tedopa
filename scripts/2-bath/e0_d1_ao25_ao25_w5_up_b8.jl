using DrWatson
@quickactivate :flucn_tedopa

ITensors.disable_warn_order()


##########################################################################
##########################################################################
function create_chain(N_chain, nb, S_pos_r, S_pos_t, n1_bsn_dim)

	tot_chain = nb * N_chain + 2  # Total number of chain sites

	b1_real_pos = (S_pos_t+1:nb:tot_chain)
	b2_real_pos = (S_pos_t+2:nb:tot_chain)
	b1_tilde_pos = (S_pos_t+3:nb:tot_chain)

	b1_real_dim = Int.(round.(range(n1_bsn_dim, stop = 2, length = length(b1_real_pos))))
	b2_real_dim = Int.(round.(range(n1_bsn_dim, stop = 2, length = length(b2_real_pos))))
	b1_tilde_dim = Int.(round.(range(n1_bsn_dim, stop = 2, length = length(b1_tilde_pos))))
	boson_dim = [0]
	for i in 1:lastindex(b1_real_dim)
		append!(boson_dim, b1_real_dim[i], b2_real_dim[i], b1_tilde_dim[i])
	end

	s_total = [(n == S_pos_r) || (n == S_pos_t) ? Index(2, "S=1/2") : Index(boson_dim[n], "Qudit") for n in 1:tot_chain]

	return tot_chain, boson_dim, b1_real_pos, b2_real_pos, b1_tilde_pos, s_total
end
##########################################################################

function create_recurr_tedopa(N_chain, nb, types, T_list, α_list, ω_C_list)

	support_cutoff = 700
	supp = (0, support_cutoff)  # Support of the weight function
	Nquad = 10^7  # Reduced number of quadrature points for speed
	N_coeff = N_chain + 1
	N_rec = 92

	ab = Vector{Any}(undef, nb)
	c_0 = Vector{Any}(undef, nb)
	for i in 1:nb
		ab[i] = Matrix{Float64}(undef, N_coeff, 2)

		T = T_list[i]
		β = 1 / T

		if types[i] == "real"
			w_fn = k -> w_fn_real(k, β, α_list[i], ω_C_list[i])
		elseif types[i] == "tilde"
			w_fn = k -> w_fn_tilde(k, β, α_list[i], ω_C_list[i])
		end

		# Calculate recurrence coefficients
		η0 = quadgk(w_fn, 0, support_cutoff)
		c_0[i] = sqrt(Complex(η0[1]))

		if N_chain >= N_rec
			ab[i][1:N_rec, 1:2] = recur_coeff(w_fn1, supp, N_rec, Nquad)
			a_100, b_100 = ab[i][N_rec, 1], ab[i][N_rec, 2]
			ab[i][N_rec+1:N_coeff, 1:2] .= repeat([a_100 b_100], N_coeff - N_rec, 1)
		else
			ab[i] .= recur_coeff(w_fn1, supp, N_coeff, Nquad)
		end
	end

	return ab[1], ab[2], ab[3], c_0[1], c_0[2], c_0[3]
end
##########################################################################

function HB(ab_list, S_pos_r, S_pos_t, types, b_pos_list, s_list)

	H = OpSum()
	for bn in 1:nb

		if types[bn] == "real"
			ω_n = ab_list[bn][1:end-1, 1]
			t_n = sqrt.(ab_list[bn][2:end, 2])
		elseif types[bn] == "tilde"
			ω_n = -1 * ab_list[bn][1:end-1, 1]
			t_n = -1 * sqrt.(ab_list[bn][2:end, 2])
		end

		for i in 1:length(b_pos_list[bn])
			j = b_pos_list[bn][i]
			H .+= (ω_n[i], "N", j)
			H .+= (t_n[i], "Adag", j, "A", j + nb)
			H .+= (t_n[i], "A", j, "Adag", j + nb)
		end
	end

	return MPO(H, s_list)
end


##########################################################################
function dw_ham(ω_0, Ω, c_0_list, ab_list, tau, S_pos, b_pos_list, s_total)
	tot_chain = length(s_total)
	S_pos_t = Int(tot_chain / 2)
	S_pos_r = S_pos_t + 1

	ω_n_REAL = ab1[1:S_pos_t-1, 1]
	ω_n_TILD = ab2[1:S_pos_t-1, 1]
	ω_n_total = append!(reverse(ω_n_TILD), [0, 0], ω_n_REAL)
	t_n_REAL = sqrt.(ab1[2:S_pos_t, 2])
	t_n_TILD = sqrt.(ab2[2:S_pos_t, 2])
	t_n_total = append!(reverse(t_n_TILD), [0, 0], t_n_REAL)

	ham = OpSum()

	for j in 2:tot_chain-1

		if j < S_pos_t
			ω_n = ω_n_total[j]
			t_n = t_n_total[j]
			ham .+= (-ω_n), "N", j
			ham .+= (-t_n), "Adag", j, "A", j - 1
			ham .+= (-t_n), "A", j, "Adag", j - 1

		elseif j == S_pos_t
			ham .-= ω_0, "Sz", j
			ham .-= Ω, "Sx", j
			ham .-= c_01, "Sx", j, "A", j - 1
			ham .-= c_01, "Sx", j, "Adag", j - 1
			ham .+= c_02, "Sx", j, "A", S_pos_r + 1
			ham .+= c_02, "Sx", j, "Adag", S_pos_r + 1

		elseif j == S_pos_r
			ham .+= ω_0, "Sz", j
			ham .+= Ω, "Sx", j
			ham .+= c_01, "Sx", j, "A", j + 1
			ham .+= c_01, "Sx", j, "Adag", j + 1
			ham .-= c_02, "Sx", j, "A", S_pos_t - 1
			ham .-= c_02, "Sx", j, "Adag", S_pos_t - 1

		else
			ω_n = ω_n_total[j]
			t_n = t_n_total[j]
			ham .+= ω_n, "N", j
			ham .+= t_n, "Adag", j, "A", j + 1
			ham .+= t_n, "A", j, "Adag", j + 1
		end
	end

	return MPO(ham, s_total)
end


##########################################################################


# Method definitions must be at the top level, not inside functions
ITensors.op(::OpName"ρ", ::SiteType"S=1/2") = [1.0 1.0; 1.0 1.0] ./ 2.0  # Adjusted normalization of the spin state
ITensors.state(::StateName"+", ::SiteType"S=1/2") = (1 / sqrt(2)) * [1; 1]  # Density matrix for qudit
ITensors.state(::StateName"up", ::SiteType"S=1/2") = [1; 0]  # Density matrix for qudit

ITensors.op(::OpName"0", ::SiteType"Qudit", d::Int) = 1.0I[1:d, 1] * 1.0I[1:d, 1]'
ITensors.state(::StateName"0", ::SiteType"Qudit", d::Int) = 1.0I[1:d, 1]
ITensors.op(::OpName"Idd", ::SiteType"Qudit", d::Int) = (1 / d) * Matrix(1.0I, d, d)
ITensors.op(::OpName"Idd", ::SiteType"S=1/2") = (1 / 2) * Matrix(1.0I, 2, 2)

let
	################# Parameters ########################################################################
	# Define filename for output
	file_name_txt_m = string(split(split(@__FILE__, ".")[end-1], string('\\'))[end], "_mQ.txt")
	file_name_txt_v = string(split(split(@__FILE__, ".")[end-1], string('\\'))[end], "_vQ.txt")
	@show file_name_txt_m
	@show file_name_txt_v

	# Define parameters for simulation
	cut = -13  # Cutoff for singular values
	cutoff = 10.0^cut
	maxdim = 50
	tau = 0.002  # Time step duration
	jump = 10  # Number of time steps for each evolution
	nt = 2500  # Number of time steps
	ttotal = nt * tau  # Total time evolution

	t_list = collect(0:1:nt) * tau


	S_pos_r = 1  # Position of the spin site
	S_pos_t = S_pos_r + 1
	N_chain = 80  # Number of chain sites for a single chain-transformed environment
	nb = 3 # Number of baths (counting real-tilde)
	types = ["real", "real", "tilde"]  # Type of baths
	n1_bsn_dim = 8  # Dimension of chain sites
	tot_chain, boson_dim, b1_real_pos, b2_real_pos, b1_tilde_pos, s_total = create_chain(N_chain, nb, S_pos_r, S_pos_t, n1_bsn_dim)

	b_pos_list = [b1_real_pos, b2_real_pos, b1_tilde_pos]
	println("Total chain length: $tot_chain")

	state = [(n == S_pos) ? "up" : "0" for n in 1:tot_chain]
	ψ = MPS(s_total, state)

	# Bath parameters
	ω_C = 5  # Bath cutoff
	ω_0 = 1  # Spin splitting
	Ω = 0
	model = (ω_0 == 1) ? "local" : "tunnel"
	model = (Ω == 1) ? "tunnel" : "local"

	T1 = 1.0  # Temperature of bath-1 real
	T2 = 1e-30  # Temperature of bath-2 real
	T_list = [T1, T2, T1]

	α = 0.25
	α_list = [α, α, α]
	ω_C_list = [ω_C, ω_C, ω_C]

	# Define functions for the weight functions
	n(ω, β) = 1 / (exp(β * ω) - 1)
	w_fn_real(k, β, α, ω_C) = (2 * α * k * exp(-k / ω_C)) * (1 + n(k, β))
	w_fn_tilde(k, β, α, ω_C) = (2 * α * k * exp(-k / ω_C)) * n(k, β)

	ab1, ab2, ab3, c_01, c_02, c_03 = create_recurr_tedopa(N_chain, nb, types, T_list, α_list, ω_C_list)

	ab_list = [ab1, ab2, ab3]
	c_0_list = [c_01, c_02, c_03]

	heat_op = HB(ab_list, S_pos, b_pos_list, s_total)
	heat_op_2 = apply(heat_op, heat_op)
	#evol = apply(dw_unit_gates(ω_0, Ω, c_01, c_02, ab1, ab2, tau, s_total), MPO(s_total, "Id"); cutoff = 1e-15)
	evol = dw_unit_gates(ω_0, Ω, c_0_list, ab_list, tau, s_total)
	ω_0, Ω, c_0_list, ab_list, tau, S_pos, b_pos_list, s_total

	# Initialize characteristic function vector
	char_fn = Vector{ComplexF64}()
	t_plot = Vector{Float64}()

	mean_Q = Float64[]
	var_Q = Float64[]

	# Time evolution of state
	U_ψ = ψ
	@show mQ = real(inner(U_ψ', heat_op, U_ψ))
	@show vQ = real(inner(heat_op, U_ψ, heat_op, U_ψ)) - mQ^2
	push!(mean_Q, mQ)
	push!(var_Q, vQ)
	write_for_loop(file_name_txt_m, string(1), "$(model) boson: T = $T, alpha = $α, N_chain = $N_chain, maxdim = $maxdim, cutoff = $cut, tau = $tau, jump = $jump, boson_dim = $n1_bsn_dim, omega = $ω_C")
	write_for_loop(file_name_txt_v, string(1), "$(model) boson: T = $T, alpha = $α, N_chain = $N_chain, maxdim = $maxdim, cutoff = $cut, tau = $tau, jump = $jump, boson_dim = $n1_bsn_dim, omega = $ω_C")
	write_for_loop(file_name_txt_m, string(2), string(mQ))
	write_for_loop(file_name_txt_v, string(2), string(vQ))

	for t in 1:nt

		U_ψ = apply(evol, U_ψ; cutoff)

		if t % jump == 0
			normalize!(U_ψ)
			ITensors.truncate!(U_ψ; maxdim = 230)
			@show mQ = real(inner(U_ψ', heat_op, U_ψ))
			@show vQ = real(inner(heat_op, U_ψ, heat_op, U_ψ)) - mQ^2
			write_for_loop(file_name_txt_m, string(t + 1), string(mQ))
			write_for_loop(file_name_txt_v, string(t + 1), string(vQ))
			push!(mean_Q, mQ)
			push!(var_Q, vQ)

			@show t * tau
			@show maxlinkdim(U_ψ)
		end

	end
	write_to_file(file_name_txt_m, "$(model) boson: T = $T, alpha = $α, N_chain = $N_chain, maxdim = $maxdim, cutoff = $cut, tau = $tau, jump = $jump, boson_dim = $n1_bsn_dim, omega = $ω_C", string(mean_Q))
	write_to_file(file_name_txt_v, "$(model) boson: T = $T, alpha = $α, N_chain = $N_chain, maxdim = $maxdim, cutoff = $cut, tau = $tau, jump = $jump, boson_dim = $n1_bsn_dim, omega = $ω_C", string(var_Q))
end
