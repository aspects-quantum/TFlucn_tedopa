using DrWatson
@quickactivate :flucn_tedopa
## https://juliadynamics.github.io/DrWatson.jl/stable/real_world/#Making-your-project-a-usable-module-1

#@disable_warn_order
ITensors.disable_warn_order()

#################            ########################################################################
################# Parameters ########################################################################
#################            ########################################################################

file_name = string(split(split(@__FILE__, ".")[end-1], string(\))[end],".txt")

ρ = [1 0;0 0]    
ρ = ρ/tr(ρ)             ## initial spin state

ITensors.op(::OpName"ρ", ::SiteType"S=1/2") = ρ


cutoff = 1E-5
#maxdim_ops = 5
maxdim = 11
tau = 10^-2             ## time step duration
nt = 2
Nbeta = 1
ttotal = nt * tau           ## TOTAL TIME evolution

N_chain = 3            ## Number of chain sites for single chain-transformed environment
tot_chain = 2*N_chain+2


# Make an array of 'site' INDICES for the (spin+chain)
S_pos_t = N_chain+1
S_pos_r = N_chain+2

println(S_pos_r)

n1_bsn_dim = 2;
b_dim = [n1_bsn_dim-round(Int64,(n1_bsn_dim-2.6)*(i-1)/(N_chain-1)) for i = 1:N_chain]           ## Dimension of chain sites
boson_dim = append!(reverse(b_dim),[0 0],b_dim)


s_total = [(n == S_pos_t) | (n == S_pos_r) ? Index(2, "S=1/2") : Index(boson_dim[n], "Qudit") for n = 1:tot_chain]
s_tot_comp = s_total[S_pos_r:tot_chain]# [(n == 1) ? Index(2, "S=1/2") : Index(boson_dim[S_pos_t + n], "Qudit") for n = 1:S_pos_t]

ITensors.op(::OpName"Idd", ::SiteType"Qudit", d::Int) = (1 / d) * Matrix(1.0I, d, d)
ITensors.op(::OpName"Idd", ::SiteType"S=1/2") = (1 / 2) * Matrix(1.0I, 2, 2)


ω_C = 5                 ## bath cutoff
ω_0 = 1                 ## spin splitting
Ω = 0                   ## independent model if Ω = 0

T = 5;                  ## temperature of bath
β = 1/T 
###

u = 0.01 # counting field parameter


## 

p=plot()

t_list = collect(0:1:nt)*tau
#α_list = [0.1,1.5]
α = .2

#cat(a, b) = reshape(append!(vec(a), vec(b)), size(a)[1:end-1]..., :)
support_cutoff = 7*10^2
supp = (0, support_cutoff)                            ## support of the weight function
Nquad = 10^7                                          ## Number of quadrature points
N_coeff = N_chain + 1

ab1 = Matrix{Float64}(undef,N_coeff,2)
ab2 = Matrix{Float64}(undef,N_coeff,2)

n(k) = 1/(exp(β*k) - 1)

#for α = α_list
w_fn1(k) = (2*α*k*exp(-k/ω_C))#*(1 + n(k))           ## weight function for real space bath
w_fn2(k) = (2*α*k*exp(-k/ω_C))#*n(k)           ## weight function for tilde space bath
η01 = quadgk(w_fn1, 0, support_cutoff)
c_01 = sqrt(η01[1])
η02 = quadgk(w_fn2, 0, support_cutoff)
c_02 = sqrt(η02[1])

if N_chain >= 90
    ab1[1:90,1:2] = recur_coeff(w_fn1, supp, 90, Nquad)    ## recurrence coefficients for for real space bath
    ab2[1:90,1:2] = recur_coeff(w_fn2, supp, 90, Nquad)    ## recurrence coefficients for for tilde space bath
    a1_100 = ab1[90,1]; b1_100 = ab1[90,2]; a2_100 = ab2[90,1]; b2_100 = ab2[90,2]
    [ab1[i,1] = a1_100  for i = 91:N_coeff]; [ab1[i,2] = b1_100  for i = 91:N_coeff]; 
    [ab2[i,1] = a2_100  for i = 91:N_coeff]; [ab2[i,2] = b2_100  for i = 91:N_coeff]; 
else
    ab1 = recur_coeff(w_fn1, supp, N_coeff, Nquad)    ## recurrence coefficients for real space bath
    ab2 = recur_coeff(w_fn2, supp, N_coeff, Nquad)    ## recurrence coefficients for tilde space bath
end



I0 = MPO(s_tot_comp, "Id")
ρ_vec = convert(MPS, I0)
for i = 1:Int(length(s_tot_comp))
    ρ_vec[i] = ρ_vec[i]*delta(s_tot_comp[i]', s_total[S_pos_t+1-i])
end
normalize!(ρ_vec)

for i = 1:Nbeta           ## perform IMAGINARY TIME evolution
    global ρ_vec = apply(exp_xHB_comp(ab1, -0.5*β/Nbeta, s_tot_comp), ρ_vec; cutoff, maxdim)
    normalize!(ρ_vec)
    println("thermal state", i)
end

@show ρ_vec 


#= 
hbb = HBB_tf(ab1, ab2, s_total)
bath_ham = MPO(s_tot_comp)
for i = 1:Int(length(s_tot_comp))
    bath_ham[S_pos_t+1-i] = hbb[i]*hbb[tot_chain+1-i]
end
@show bath_ham

htot = HTOT_tf(ω_0, Ω, c_01, c_02, ab1, ab2, s_total)
tot_ham = MPO(s_tot_comp)
for i = 1:Int(length(s_tot_comp))
    tot_ham[S_pos_t+1-i] = htot[i]*htot[tot_chain+1-i]
end
 =#

count_gates = exp_xHB_comp(ab1, -im*u, s_tot_comp)
evol = tot_gate(ω_0, Ω, c_01, c_02, ab1, ab2, tau, s_total)

char_fn=Vector{ComplexF64}()

cg_ρ = apply(count_gates, ρ_vec; cutoff, maxdim)
U_cg_ρ = cg_ρ
U_ρ = ρ_vec
cg_U_ρ = cg_ρ

chi = inner(cg_U_ρ, U_cg_ρ)
println(chi[1])
push!(char_fn,chi[1])

for t in 1:nt
    global U_cg_ρ = apply(evol, U_cg_ρ; cutoff, maxdim)
    global U_ρ = apply(evol, U_ρ; cutoff, maxdim)
    global cg_U_ρ = apply(count_gates, U_ρ; cutoff, maxdim)

    chi = inner(cg_U_ρ, U_cg_ρ)
    println(chi[1])
    write_for_loop(file_name, string(t), string(chi[1]))
    push!(char_fn,chi[1])
end
chi_pu = char_fn

mean_Q = imag.(chi_pu)/u
@show real(mean_Q)
write_to_file(file_name, string("Independent boson: T = 5, alpha = .2, N_chain = 30, maxdim = 110, tau = 0.01"), string(real(mean_Q)))

#plot!(t_list,real(mean_Q),label= "α = $α")

plot!((collect(0:1:length(mean_Q)-1)*tau),real(mean_Q),label= "α = $α")
#end


plot!(legend=:topright)
xlabel!("t")
title = string("<Q>, N_ch = ", N_chain,", b_dim = ", boson_dim,", u = ",u, ", k_max = ", support_cutoff)
title!(title)
display("image/png", p)


a = "tebd_ind_T5_ao2_ho01.png";
#savefig(a)
safesave(a,p)
   
#= 
using HDF5
f = h5open("myfile.h5","r")
T = read(f,"thermal_state",MPS)
close(f) =#