using DrWatson
using Base.Threads

@quickactivate :flucn_tedopa

BLAS.set_num_threads(1)
ITensors.disable_warn_order()


##########################################################################
##########################################################################
function create_chain(N_chain, nb, S_pos_r, S_pos_t, n1_bsn_dim)

	tot_chain = nb * N_chain + 2  # Total number of chain sites

	b1_real_pos = (S_pos_r+1:nb:tot_chain)
	b2_real_pos = (S_pos_r+2:nb:tot_chain)
	b1_tilde_pos = (S_pos_r+3:nb:tot_chain)

	b1_real_dim = Int.(round.(range(n1_bsn_dim[1], stop = 2, length = length(b1_real_pos))))
	b2_real_dim = Int.(round.(range(n1_bsn_dim[2], stop = 2, length = length(b2_real_pos))))
	b1_tilde_dim = Int.(round.(range(n1_bsn_dim[3], stop = 2, length = length(b1_tilde_pos))))
	boson_dim = [0, 0]
	for i in 1:lastindex(b1_real_dim)
		append!(boson_dim, b1_real_dim[i], b2_real_dim[i], b1_tilde_dim[i])
	end

	s_total = [(n == S_pos_r) || (n == S_pos_t) ? Index(2, "S=1/2") : Index(boson_dim[n], "Qudit") for n in 1:tot_chain]

	return tot_chain, boson_dim, b1_real_pos, b2_real_pos, b1_tilde_pos, s_total
end
##########################################################################

function create_recurr_tedopa(N_chain, nb, types, T_list, α_list, ω_C_list)

	# Define functions for the weight functions
	n(ω, β) = 1 / (exp(β * ω) - 1)
	w_fn_real(k, β, α, ω_C) = (2 * α * k * exp(-k / ω_C)) * (1 + n(k, β))
	w_fn_tilde(k, β, α, ω_C) = (2 * α * k * exp(-k / ω_C)) * n(k, β)

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
		@show c_0[i] = sqrt(Complex(η0[1]))

		if N_chain >= N_rec
			ab[i][1:N_rec, 1:2] = recur_coeff(w_fn, supp, N_rec, Nquad)
			a_100, b_100 = ab[i][N_rec, 1], ab[i][N_rec, 2]
			ab[i][N_rec+1:N_coeff, 1:2] .= repeat([a_100 b_100], N_coeff - N_rec, 1)
		else
			ab[i] .= recur_coeff(w_fn, supp, N_coeff, Nquad)
		end
	end

	return ab[1], ab[2], ab[3], c_0[1], c_0[2], c_0[3]
end
##########################################################################

function HB(which_baths, ab_list, types, b_pos_list, s_list)

	H = OpSum()
	for bn in which_baths

		if types[bn] == "real"
			ω_n = ab_list[bn][1:end-1, 1]
			t_n = sqrt.(ab_list[bn][2:end, 2])
		elseif types[bn] == "tilde"
			ω_n = -1 * ab_list[bn][1:end-1, 1]
			t_n = -1 * sqrt.(ab_list[bn][2:end, 2])
		end

		for i in 1:lastindex(b_pos_list[bn])-1
			j = b_pos_list[bn][i]
			k = b_pos_list[bn][i+1]
			H .+= (ω_n[i], "N", j)
			H .+= (t_n[i], "Adag", j, "A", k)
			H .+= (t_n[i], "A", j, "Adag", k)
		end
	end

	return MPO(H, s_list)
end


##########################################################################
function dw_ham(ω_0, Ω, c_0_list, L, ab_list, S_pos_r, S_pos_t, nb, types, b_pos_list, s_list)
	ham = OpSum()

	ham .+= (ω_0, "Sz", S_pos_r)
	ham .+= (Ω, "Sx", S_pos_r)
	for bn in 1:lastindex(c_0_list)
		c_0 = c_0_list[bn]
		if c_0 == 0
			continue
		else
			j = b_pos_list[bn][1]
			if types[bn] == "real"
				ham .+= (c_0, L, S_pos_r, "A", j)
				ham .+= (c_0, L, S_pos_r, "Adag", j)
			elseif types[bn] == "tilde"
				ham .+= (c_0, L, S_pos_r, "A", j)
				ham .+= (c_0, L, S_pos_r, "Adag", j)
			end
		end
	end

	c0_uvrev_list = circshift(c_0_list, -Int(lastindex(c_0_list) / 2))
	ham .+= (-ω_0, "Sz", S_pos_t)
	ham .+= (-Ω, "Sx", S_pos_t)
	bn = 1
	for m in 1:lastindex(c_0_list)
		c_0 = c0_uvrev_list[m]
		if c_0 == 0
			continue
		else
			j = b_pos_list[bn][1]
			if types[bn] == "real"
				ham .+= (-c_0, L, S_pos_t, "Adag", j)
				ham .+= (-c_0, L, S_pos_t, "A", j)
			elseif types[bn] == "tilde"
				ham .+= (-c_0, L, S_pos_t, "Adag", j)
				ham .+= (-c_0, L, S_pos_t, "A", j)
			end
			bn += 1
		end
	end

	for bn in 1:nb

		if types[bn] == "real"
			ω_n = ab_list[bn][1:end-1, 1]
			t_n = sqrt.(ab_list[bn][2:end, 2])
		elseif types[bn] == "tilde"
			ω_n = -1 * ab_list[bn][1:end-1, 1]
			t_n = -1 * sqrt.(ab_list[bn][2:end, 2])
		end
		for i in 1:length(b_pos_list[bn])-1
			j = b_pos_list[bn][i]
			k = b_pos_list[bn][i+1]
			ham .+= (ω_n[i], "N", j)
			ham .+= (t_n[i], "Adag", j, "A", k)
			ham .+= (t_n[i], "A", j, "Adag", k)
		end
	end

	return MPO(ham, s_list; splitblocks = true)
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
	file_name_txt_m = string(split(split(@__FILE__, ".")[end-1], string('\\'))[end], "_mJ.txt")
	file_name_txt_v = string(split(split(@__FILE__, ".")[end-1], string('\\'))[end], "_vJ.txt")
	file_name_txt_N = string(split(split(@__FILE__, ".")[end-1], string('\\'))[end], "_N.txt")

	# Define parameters for simulation
	cut = -16 # Cutoff for singular values
	cutoff = 10.0^cut
	maxdim = 100
	tau = 0.001  # Time step duration
	jump = 10 # Number of time steps between recorded data
	nt = 1000  # Number of time steps
	ttotal = nt * tau  # Total time evolution
	tdvp_steps = 1 # Number of substeps in each tdvp step

	###############################################################################################
	# Create system
	S_pos_r = 2  # Position of the spin site
	S_pos_t = S_pos_r - 1
	N_chain = 100  # Number of chain sites for a single chain-transformed environment
	nb = 3 # Number of baths (counting real-tilde)
	types = ["real", "real", "tilde", "tilde"]  # Type of baths
	#n1_bsn_dim = 8  # Dimension of chain sites
	n1_bsn_dim = [15, 15, 12]
	tot_chain, boson_dim, b1_real_pos, b2_real_pos, b1_tilde_pos, s_total = create_chain(N_chain, nb, S_pos_r, S_pos_t, n1_bsn_dim)

	b_pos_list = [b1_real_pos, b2_real_pos, b1_tilde_pos]

	state_base = [(n == S_pos_r) || (n == S_pos_t) ? "up" : "0" for n in 1:tot_chain]
	ψ_base = MPS(s_total, state_base)

	L = "Sx"

	# Bath parameters
	ω_C = 0.25  # Bath cutoff
	ω_0 = 1  # Spin splitting
	Ω = 0
	model = (ω_0 == 1) ? "local" : "tunnel"
	model = (Ω == 1) ? "tunnel" : "local"

	T1 = 1.0  # Temperature of bath-1 real
	T2 = 1e-30  # Temperature of bath-2 real
	T_list = [T1, T2, T1, T2]

	α1 = 0.25
	α2 = 1.25
	α_list = [α1, α2, α1, α2]
	ω_C_list = [ω_C, ω_C, ω_C, ω_C]

	ab1, ab2, ab3, c_01, c_02, c_03 = create_recurr_tedopa(N_chain, nb, types, T_list, α_list, ω_C_list)

	ab_list = [ab1, ab2, ab3]
	c_0_list = [c_01, c_02, c_03, 0]


	mean_J = Float32[]
	var_J = Float32[]

	N_temp = 30  # Temporary chain length 
	println("Total chain length: $(2+nb*N_temp)")

	s_list = s_total[1:2+nb*N_temp]
	b_pos_temp = [b1_real_pos[1:N_temp], b2_real_pos[1:N_temp], b1_tilde_pos[1:N_temp]]

	state = [(n == S_pos_r) || (n == S_pos_t) ? "up" : "0" for n in 1:2+3N_temp] # Initial state
	ψ = MPS(s_list, state)


	which_baths = [1, 3]
	heat_op1 = HB(which_baths, ab_list, types, b_pos_temp, s_list)
	which_baths = [2]
	heat_op2 = HB(which_baths, ab_list, types, b_pos_temp, s_list)

	J0 = heat_op2 - heat_op1
	J0_dag = noprime(linkinds, swapprime(dag(J0), 0 => 1))
	J = 0.5 * add(J0, J0_dag; cutoff = 1e-17)

	evol = dw_ham(ω_0, Ω, c_0_list, L, ab_list, S_pos_r, S_pos_t, nb, types, b_pos_temp, s_list)
	#= @show maxlinkdim(evol0)
	evol0_dag = noprime(linkinds, swapprime(dag(evol0), 0 => 1))
	evol = 0.5 * add(evol0, evol0_dag; cutoff = 1e-17) =#
	@show maxlinkdim(evol)
	orthogonalize!(evol, S_pos_r)

	# Time evolution of state
	U_ψ = ψ
	orthogonalize!(U_ψ, S_pos_r)

	mJ = real(inner(U_ψ', J, U_ψ)) / tau
	vJ = real(inner(J, U_ψ, J, U_ψ)) - mJ^2 / tau
	push!(mean_J, mJ)
	push!(var_J, vJ)
	#= @show mQ2 = real(inner(U_ψ', heat_op2, U_ψ))
	@show vQ2 = real(inner(heat_op2, U_ψ, heat_op2, U_ψ)) - mQ2^2
	push!(mean_Q2, mQ2)
	push!(var_Q2, vQ2) =#


	maxdim_1site = 200
	nsites = 2


	write_for_loop(file_name_txt_m, string(1), "$(model) boson: T = [$T1, $T2], alpha = [$α1, $α2], N_chain = $N_temp, maxdim = $maxdim_1site, cutoff = $cut, tau = $tau, jump = $jump, boson_dim = $n1_bsn_dim, omega = $ω_C")
	write_for_loop(file_name_txt_v, string(1), "$(model) boson: T = [$T1, $T2], alpha = [$α1, $α2], N_chain = $N_temp, maxdim = $maxdim_1site, cutoff = $cut, tau = $tau, jump = $jump, boson_dim = $n1_bsn_dim, omega = $ω_C")
	write_for_loop(file_name_txt_N, string(1), "$(model) boson: T = [$T1, $T2], alpha = [$α1, $α2], N_chain = $N_temp, maxdim = $maxdim_1site, cutoff = $cut, tau = $tau, jump = $jump, boson_dim = $n1_bsn_dim, omega = $ω_C")
	write_for_loop(file_name_txt_m, string(2), string(mJ))
	write_for_loop(file_name_txt_v, string(2), string(vJ))

	# coarse_future will be created *after* the first checkpoint (pipelined)
	coarse_future = nothing
	ψ_chk = deepcopy(U_ψ)  # start-of-interval checkpoint
	last_t = 0
	@assert iseven(jump) "jump must be even for coarse 2*dt"

	for t in 1:nt
		@show t * tau
		@show maxlinkdim(U_ψ)

		# fine (live) evolution with dt = tau
		U_ψ = tdvp(evol, -1im * tau, U_ψ; nsteps = tdvp_steps, nsite = nsites,
			normalize = true, cutoff = cutoff, maxdim = maxdim_1site)
		# orthogonalize occasionally; doing it every step slows things down
		orthogonalize!(U_ψ, S_pos_r)

		if t % jump == 0

			if (iseven(nsites) == true) || (maxlinkdim(U_ψ) >= maxdim_1site)
				nsites = 1
			else
				nsites = 2
			end

			# ---- Fine (A) from live state ----
			J_U_ψ = apply(J, U_ψ; cutoff = 1e-17)
			e1A   = real(inner(U_ψ, J_U_ψ))
			e2A   = real(inner(J_U_ψ, J_U_ψ))  # = ||Jψ||^2
			mA    = e1A / t
			vA    = (e2A - e1A^2) / t

			# ---- Coarse (B) for this interval ----
			Δsteps = t - last_t
			@assert iseven(Δsteps)
			if coarse_future === nothing
				# First interval: run coarse synchronously (no pre-launch stall)
				ψB = deepcopy(ψ_chk)
				@inbounds for _ in 1:(Δsteps÷2)
					ψB = tdvp(evol, -1im * (2 * tau), ψB; nsteps = tdvp_steps, nsite = 2,
						normalize = true, cutoff = cutoff, maxdim = maxdim_1site)
				end
				JψB = apply(J, ψB; cutoff = 1e-17)
				e1B = real(inner(ψB, JψB))
				e2B = real(inner(JψB, JψB))
				mB  = e1B / t
				vB  = (e2B - e1B^2) / t
			else
				mB, vB = fetch(coarse_future)
			end

			# ---- Richardson (p=1): 2A - B ----
			@show mJ = (2 * mA - mB) / tau
			@show vJ = (2 * vA - vB) / tau
			push!(mean_J, mJ)
			push!(var_J, vJ)
			write_for_loop(file_name_txt_m, string(t + 1), string(mJ))
			write_for_loop(file_name_txt_v, string(t + 1), string(vJ))

			# ---- Pipeline: update checkpoint and spawn the NEXT coarse run ----
			ψ_chk  = deepcopy(U_ψ)
			last_t = t
			if t + jump <= nt
				local Tsteps_target = t + jump
				local ψ_start = ψ_chk
				coarse_future = @spawn begin
					ψB2 = deepcopy(ψ_start)
					@inbounds for _ in 1:(jump÷2)
						ψB2 = tdvp(evol, -1im * (2 * tau), ψB2; nsteps = tdvp_steps, nsite = nsites,
							normalize = true, cutoff = cutoff, maxdim = maxdim_1site)
					end
					JψB2 = apply(J, ψB2; cutoff = 1e-17)
					e1B2 = real(inner(ψB2, JψB2))
					e2B2 = real(inner(JψB2, JψB2))
					(e1B2 / Tsteps_target, (e2B2 - e1B2^2) / Tsteps_target)
				end
			end
			#= n_list = []
			for i in 1:nb*N_temp
				ni = MPO(OpSum() + (1, "N", 2 + i), s_list)
				push!(n_list, real(inner(U_ψ', ni, U_ψ)))
			end =#
			χ_list = linkdims(U_ψ)
			@show i = findlast(>(1), χ_list)
			isnothing(i) ? lastindex(s_list) : χ_list[i]
			write_for_loop(file_name_txt_N, string(t + 1), string(χ_list))
		end
	end

	write_to_file(file_name_txt_m, "$(model) boson: T = [$T1, $T2], alpha = [$α1, $α2], N_chain = $N_temp, maxdim = $maxdim_1site, cutoff = $cut, tau = $tau, jump = $jump, boson_dim = $n1_bsn_dim, omega = $ω_C", string(mean_J))
	write_to_file(file_name_txt_v, "$(model) boson: T = [$T1, $T2], alpha = [$α1, $α2], N_chain = $N_temp, maxdim = $maxdim_1site, cutoff = $cut, tau = $tau, jump = $jump, boson_dim = $n1_bsn_dim, omega = $ω_C", string(var_J))
end
