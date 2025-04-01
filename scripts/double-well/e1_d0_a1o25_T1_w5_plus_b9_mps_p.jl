using DrWatson
@quickactivate :flucn_tedopa

ITensors.disable_warn_order()


##########################################################################
##########################################################################

function HB(ab1, ab2, s_total)

	tot_chain = length(s_total)
	S_pos = Int((tot_chain + 1) / 2)

	ω_n_REAL = ab1[1:S_pos-1, 1]
	ω_n_TILD = ab2[1:S_pos-1, 1]
	ω_n_total = append!(reverse(ω_n_TILD), [0], ω_n_REAL)
	t_n_REAL = sqrt.(ab1[2:S_pos, 2])
	t_n_TILD = sqrt.(ab2[2:S_pos, 2])
	t_n_total = append!(reverse(t_n_TILD), [0], t_n_REAL)

	H = OpSum()

	for j in 2:tot_chain-1

		if j < S_pos
			ω_n = ω_n_total[j]
			t_n = t_n_total[j]

			H .+= (-ω_n, "N", j, "Id", j - 1)
			H .+= (-t_n, "Adag", j, "A", j - 1)
			H .+= (-t_n, "A", j, "Adag", j - 1)
		elseif j == S_pos
			continue

		else
			ω_n = ω_n_total[j]
			t_n = t_n_total[j]

			H .+= (ω_n, "N", j, "Id", j + 1)
			H .+= (t_n, "Adag", j, "A", j + 1)
			H .+= (t_n, "A", j, "Adag", j + 1)
		end
	end

	return MPO(H, s_total)
end
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
				 c_01 * op("Sx", s_total[j]) * op("Adag", s_total[j+1])
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
	gates_rev = reverse(gates)
	return append!(gates, gates_rev)

end


##########################################################################


# Method definitions must be at the top level, not inside functions
ITensors.op(::OpName"ρ", ::SiteType"S=1/2") = [1.0 1.0; 1.0 1.0] ./ 2.0  # Adjusted normalization of the spin state
ITensors.state(::StateName"+", ::SiteType"S=1/2") = (1 / sqrt(2)) * [1; 1]  # Density matrix for qudit
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

	N_chain = 180  # Number of chain sites for a single chain-transformed environment
	tot_chain = 2 * N_chain + 1  # Total number of chain sites
	S_pos = N_chain + 1  # Position of the spin site

	println(S_pos)

	n1_bsn_dim = 9  # Dimension of chain sites
	b_dim_real = [n1_bsn_dim - round(Int64, (n1_bsn_dim - 1.6) * (i - 1) / (N_chain - 1)) for i in 1:N_chain]  # Dimension of chain sites
	b_dim_tilde = [Int(n1_bsn_dim - 4 - round(Int64, (n1_bsn_dim - 4 - 1.6) * (i - 1) / (N_chain - 1))) for i in 1:N_chain]
	boson_dim = append!(reverse(b_dim_tilde), [0], b_dim_real)

	s_total = [(n == S_pos) ? Index(2, "S=1/2") : Index(boson_dim[n], "Qudit") for n in 1:tot_chain]
	
	state = [(n == S_pos) ? "+" : "0" for n in 1:tot_chain]
	ψ = MPS(s_total, state)

	# Bath parameters
	ω_C = 5  # Bath cutoff
	ω_0 = 1  # Spin splitting
	Ω = 0
	model = (ω_0 == 1) ? "local" : "tunnel"
	model = (Ω == 1) ? "tunnel" : "local"

	T = 1.  # Temperature of bath
	β = 1 / T

	t_list = collect(0:1:nt) * tau
	α = 1.25
	mean_Q = Float64[]

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


	heat_op = HB(ab1, ab2, s_total)
	heat_op_2 = apply(heat_op, heat_op)
	#evol = apply(dw_unit_gates(ω_0, Ω, c_01, c_02, ab1, ab2, tau, s_total), MPO(s_total, "Id"); cutoff = 1e-15)
	evol = dw_unit_gates(ω_0, Ω, c_01, c_02, ab1, ab2, tau, s_total)

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
			ITensors.truncate!(U_ψ; maxdim = 200)
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
