using DrWatson
@quickactivate :flucn_tedopa

ITensors.disable_warn_order()



# Method definitions must be at the top level, not inside functions
ITensors.op(::OpName"ρ", ::SiteType"S=1/2") = [1 -1; -1 1] / 2  # Adjusted normalization of the spin state
ITensors.op(::OpName"0", ::SiteType"Qudit", d::Int) = 1.0I[1:d, 1] * 1.0I[1:d, 1]'
ITensors.op(::OpName"Idd", ::SiteType"Qudit", d::Int) = (1 / d) * Matrix(1.0I, d, d)
ITensors.op(::OpName"Idd", ::SiteType"S=1/2") = (1 / 2) * Matrix(1.0I, 2, 2)

let
	################# Parameters ########################################################################
	# Define filename for output
	file_name_txt = string(split(split(@__FILE__, ".")[end-1], string('\\'))[end], ".txt")
	@show file_name_txt

	# Define parameters for simulation
	cut = -11
	cutoff = 10.0^cut
	maxdim = 200
	tau = 0.002  # Time step duration
	nt = 700
	ttotal = nt * tau  # Total time evolution

	N_chain = 160  # Number of chain sites for a single chain-transformed environment
	tot_chain = 2 * N_chain + 1
	S_pos = N_chain + 1

	println(S_pos)

	n1_bsn_dim = 11
	b_dim_real = [n1_bsn_dim - round(Int64, (n1_bsn_dim - 3.6) * (i - 1) / (N_chain - 1)) for i in 1:N_chain]  # Dimension of chain sites
	b_dim_tilde = [Int(n1_bsn_dim - 5 - round(Int64, (n1_bsn_dim - 5 - 1.6) * (i - 1) / (N_chain - 1))) for i in 1:N_chain]
	boson_dim = append!(reverse(b_dim_tilde), [0], b_dim_real)

	s_total = [(n == S_pos) ? Index(2, "S=1/2") : Index(boson_dim[n], "Qudit") for n in 1:tot_chain]

	state = [(n == S_pos) ? "ρ" : "0" for n in 1:tot_chain]
	ρ = normalize(MPO(s_total, state))

	# Bath parameters
	ω_C = 5  # Bath cutoff
	ω_0 = 0  # Spin splitting
	Ω = 1  # Independent model if Ω = 0
	model = (ω_0 == 1) ? "independent" : "unbiased"

	T = .1  # Temperature of bath
	β = 1 / T
	u = 0.005  # Counting field parameter

	t_list = collect(0:1:nt) * tau
	α = 0.5
	mean_Q = Float64[]

	support_cutoff = 700
	supp = (0, support_cutoff)  # Support of the weight function
	Nquad = 10^7  # Reduced number of quadrature points for speed
	N_coeff = N_chain + 1
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

	# Time evolution of density matrix
	U_ρ_Ud = copy(ρ)
	cd_U_ρ_Ud = apply(count_gates, U_ρ_Ud; cutoff = cutoff * 0.01)
	chi = tr(cd_U_ρ_Ud)
	push!(char_fn, chi[1])
	write_for_loop(file_name_txt, string(1), "$(model) boson: T = $T, alpha = $α, N_chain = $N_chain, maxdim = $maxdim, cutoff = $cut, tau = $tau, u = $u, boson_dim = $n1_bsn_dim")

	for t in 1:nt
		U_ρ_Ud = apply(evol, U_ρ_Ud; cutoff, maxdim, apply_dag = true)
		U_ρ_Ud = add(0.5 * swapprime(dag(U_ρ_Ud), 0 => 1), 0.5 * U_ρ_Ud; maxdim = 200)
		U_ρ_Ud = normalize(ITensors.truncate(ITensors.truncate(U_ρ_Ud; maxdim = 10, site_range = (1:(Int(floor(3*N_chain / 4))))); maxdim = 20, site_range = (tot_chain:tot_chain-(Int(floor((3*N_chain / 4)))))))
		if t % 10 == 0
			orthogonalize!(U_ρ_Ud, findfirst(==(maxlinkdim(U_ρ_Ud)), linkdims(U_ρ_Ud)))
			cd_U_ρ_Ud = apply(count_gates, U_ρ_Ud; cutoff = cutoff * 0.01)
			chi = tr(cd_U_ρ_Ud)
			write_for_loop(file_name_txt, string(t + 1), string(imag(chi[1]) / u))
			push!(char_fn, chi[1])
			@show imag(chi[1]) / u
		end
		@show t * tau
		@show maxlinkdim(U_ρ_Ud)
	end

	# Compute and plot mean_Q
	mean_Q = imag.(char_fn) / u
	
	write_to_file(file_name_txt, "$(model) boson: T = $T, alpha = $α, N_chain = $N_chain, maxdim = $maxdim, cutoff = $cut, tau = $tau, u = $u, boson_dim = $n1_bsn_dim", string(real(mean_Q)))
end
