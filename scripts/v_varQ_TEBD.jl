using DrWatson
@quickactivate :flucn_tedopa
## https://juliadynamics.github.io/DrWatson.jl/stable/real_world/#Making-your-project-a-usable-module-1

#@disable_warn_order
ITensors.disable_warn_order()

#################            ########################################################################
################# Parameters ########################################################################
#################            ########################################################################


file_name_txt = string(split(split(@__FILE__, ".")[end-1], string(\))[end],".txt")
file_name_png = string(split(split(@__FILE__, ".")[end-1], string(\))[end],".png")


ρ = [1 0;0 0]    
ρ = ρ/tr(ρ)             ## initial spin state

ITensors.op(::OpName"ρ", ::SiteType"S=1/2") = ρ

cut = -8
cutoff = 10.0^cut
#maxdim_ops = 5
maxdim = 150
tau = 10^-2             ## time step duration
nt = 50
Nbeta = 80
ttotal = nt * tau           ## TOTAL TIME evolution

N_chain = 15            ## Number of chain sites for single chain-transformed environment
tot_chain = 2*N_chain+2


# Make an array of 'site' INDICES for the (spin+chain)
S_pos_t = N_chain+1
S_pos_r = N_chain+2

println(S_pos_r)

n1_bsn_dim = 5;
b_dim = [n1_bsn_dim-round(Int64,(n1_bsn_dim-2.6)*(i-1)/(N_chain-1)) for i = 1:N_chain]           ## Dimension of chain sites
boson_dim = append!(reverse(b_dim),[0 0],b_dim)


s_total = [(n == S_pos_t) | (n == S_pos_r) ? Index(2, "S=1/2") : Index(boson_dim[n], "Qudit") for n = 1:tot_chain]
s_tot_comp = s_total[S_pos_r:tot_chain]# [(n == 1) ? Index(2, "S=1/2") : Index(boson_dim[S_pos_t + n], "Qudit") for n = 1:S_pos_t]

ITensors.op(::OpName"Idd", ::SiteType"Qudit", d::Int) = (1 / d) * Matrix(1.0I, d, d)
ITensors.op(::OpName"Idd", ::SiteType"S=1/2") = (1 / 2) * Matrix(1.0I, 2, 2)


ω_C = 5                 ## bath cutoff
ω_0 = 1                 ## spin splitting
Ω = 0                   ## independent model if Ω = 0

if ω_0 == 1 
    model = "independent"
else
    model = "unbiased"
end

T = 1;                  ## temperature of bath
β = 1/T 
###

u = 0.01 # counting field parameter

## 

p=plot()

t_list = collect(0:1:nt)*tau
#α_list = [0.1,1.5]
α = .1

#cat(a, b) = reshape(append!(vec(a), vec(b)), size(a)[1:end-1]..., :)
support_cutoff = 7*10^2
supp = (0, support_cutoff)                            ## support of the weight function
Nquad = 10^7                                          ## Number of quadrature points
N_coeff = N_chain + 1

ab = Matrix{Float64}(undef,N_coeff,2)

#for α = α_list
w_fn(k) = (2*α*k*exp(-k/ω_C))#*(1 + n(k))           ## weight function for real space bath
η0 = quadgk(w_fn, 0, support_cutoff)
c_0 = sqrt(η0[1])

if N_chain >= 90
    ab[1:90,1:2] = recur_coeff(w_fn, supp, 90, Nquad)    ## recurrence coefficients for for real space bath
    a1_100 = ab[90,1]; b1_100 = ab[90,2]
    [ab[i,1] = a_100  for i = 91:N_coeff]; [ab[i,2] = b_100  for i = 91:N_coeff]; 
else
    ab = recur_coeff(w_fn, supp, N_coeff, Nquad)    ## recurrence coefficients for real space bath
end


I0 = MPO(s_tot_comp, "Id")
ρ_vec = convert(MPS, I0)
for i = 1:Int(length(s_tot_comp))
    ρ_vec[i] = ρ_vec[i]*delta(s_tot_comp[i]', s_total[S_pos_t+1-i])
end
normalize!(ρ_vec)

for i = 1:Nbeta           ## perform IMAGINARY TIME evolution
    global ρ_vec = apply(exp_xHB_comp(ab, -0.5*β/Nbeta, s_total), ρ_vec; cutoff, maxdim)
    normalize!(ρ_vec)
    println("thermal state", i)
end

@show ρ_vec 


count_gates_1 = exp_xHB_comp(ab, -im*u, s_total)
count_gates_2 = exp_xHB_comp(ab, -im*2*u, s_total)
evol = tot_gate(ω_0, Ω, c_0, ab, tau, s_total)

char_fn_1=Vector{ComplexF64}()
char_fn_2=Vector{ComplexF64}()

cg_ρ_1 = apply(count_gates_1, ρ_vec; cutoff, maxdim)
cg_ρ_2 = apply(count_gates_2, ρ_vec; cutoff, maxdim)
U_cg_ρ_1 = cg_ρ_1
U_cg_ρ_2 = cg_ρ_2
U_ρ = ρ_vec
cg_U_ρ_1 = cg_ρ_1
cg_U_ρ_2 = cg_ρ_2

chi_1 = inner(cg_U_ρ_1, U_cg_ρ_1)
chi_2 = inner(cg_U_ρ_2, U_cg_ρ_2)

println(real(-(log.(chi_2[1])-2*log.(chi_1[1]))/(u^2)))
push!(char_fn_1,chi_1[1])
push!(char_fn_2,chi_2[1])


write_for_loop(file_name_txt, string(1), "Indep boson: T = $T, alpha = $α, N_chain = $N_chain, maxdim = $maxdim, cutoff = $cut, tau = $tau")
for t in 1:nt
    global U_cg_ρ_1 = apply(evol, U_cg_ρ_1; cutoff, maxdim)
    global U_ρ = apply(evol, U_ρ; cutoff, maxdim)
    global cg_U_ρ_1 = apply(count_gates_1, U_ρ; cutoff, maxdim)
    global U_cg_ρ_2 = apply(evol, U_cg_ρ_2; cutoff, maxdim)
    global U_ρ = apply(evol, U_ρ; cutoff, maxdim)
    global cg_U_ρ_2 = apply(count_gates_2, U_ρ; cutoff, maxdim)
    chi_1 = inner(cg_U_ρ_1, U_cg_ρ_1)
    chi_2 = inner(cg_U_ρ_2, U_cg_ρ_2)
    println(real(-(log.(chi_2[1])-2*log.(chi_1[1]))/(u^2)))
    push!(char_fn_1,chi_1[1])
    push!(char_fn_2,chi_2[1])
    write_for_loop(file_name_txt, string(t+1), string(real(-(log.(chi_2[1])-2*log.(chi_1[1]))/(u^2))))
end
chi_pu = char_fn_1
chi_2pu = char_fn_2

var_Q = -(log.(chi_2pu)-2*log.(chi_pu))/(u^2)
@show real(var_Q)

#plot!(t_list,real(mean_Q),label= "α = $α")
#= 
plot!((collect(0:1:length(mean_Q)-1)*tau),real(mean_Q),label= "α = $α")
#end


plot!(legend=:topright)
xlabel!("t")
title = string("<Q>, N_ch = ", N_chain,", b_dim = ", boson_dim,", u = ",u, ", k_max = ", support_cutoff)
title!(title)
display("image/png", p)


a = "tebd_ind_T5_ao2_ho01.png";
savefig(a)
    =#
#= 
using HDF5
f = h5open("myfile.h5","r")
T = read(f,"thermal_state",MPS)
close(f) =#