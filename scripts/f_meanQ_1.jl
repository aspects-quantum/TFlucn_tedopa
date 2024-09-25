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

ρ0 = [1 0;0 0]    
ρ0 = ρ0/tr(ρ0)             ## initial spin state

ITensors.op(::OpName"ρ", ::SiteType"S=1/2") = ρ0



cut = -8
cutoff = 10.0^cut
#maxdim_ops = 5
maxdim = 1000
tau = 3*10^-3             ## time step duration
nt = 500
ttotal = nt * tau           ## TOTAL TIME evolution

N_chain = 90           ## Number of chain sites for single chain-transformed environment
tot_chain = 2*N_chain+2
S_pos_t = N_chain+1
S_pos_r = N_chain+2

println(S_pos_r)

n1_bsn_dim = 6;
b_dim = [n1_bsn_dim-round(Int64,(n1_bsn_dim-2.6)*(i-1)/(N_chain-1))    for i = 1:N_chain]   #     ## Dimension of chain sites
boson_dim = append!(reverse(b_dim),[0 0],b_dim)

s_total = [(n == S_pos_t) | (n == S_pos_r) ? Index(2, "S=1/2") : Index(boson_dim[n], "Qudit") for n = 1:tot_chain]
s_tot_comp = s_total[S_pos_r:tot_chain]# [(n == 1) ? Index(2, "S=1/2") : Index(boson_dim[S_pos_t + n], "Qudit") for n = 1:S_pos_t]


ITensors.op(::OpName"0", ::SiteType"Qudit", d::Int) = 1.0I[1:d,1]*1.0I[1:d,1]'
ITensors.op(::OpName"Idd", ::SiteType"Qudit", d::Int) = (1 / d) * Matrix(1.0I, d, d)
ITensors.op(::OpName"Idd", ::SiteType"S=1/2") = (1 / 2) * Matrix(1.0I, 2, 2)


ω_C = 5                 ## bath cutoff
ω_0 = 1                 ## spin splitting
Ω = 0                   ## independent model if Ω = 0

T = 1;                  ## temperature of bath
β = 1/T 
###

u = 0.01 # counting field parameter

t_list = collect(0:1:nt)*tau
#α_list = [0.1,1.5]
α = .2

#cat(a, b) = reshape(append!(vec(a), vec(b)), size(a)[1:end-1]..., :)
support_cutoff = 7*10^3
supp = (0, support_cutoff)                            ## support of the weight function
Nquad = 10^7                                          ## Number of quadrature points
N_coeff = N_chain + 1

ab1 = Matrix{Float64}(undef,N_coeff,2)
ab2 = Matrix{Float64}(undef,N_coeff,2)

n(ω) = 1/(exp(β*ω) - 1)

#for α = α_list

w_fn1(k) = (2*α*k*exp(-k/ω_C))*(1 + n(k))         ## weight function for real space bath
w_fn2(k) = (2*α*k*exp(-k/ω_C))*(n(k))                ## weight function for tilde space bath
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
    ab2 = recur_coeff(w_fn2, supp, N_coeff, Nquad) 
end

state = [(n == S_pos_r) | (n == S_pos_t) ? "ρ" : "0" for n = 1:tot_chain]
ρ = MPO(s_total,state)

count_gates = exp_xHB_comp(ab1, ab2, im*u, s_total)
evol = tot_gate(ω_0, Ω, c_01, c_02, ab1, ab2, tau, s_total)

char_fn=Vector{ComplexF64}()

U_ρ_Ud = ρ
cd_U_ρ_Ud = apply(count_gates, U_ρ_Ud; cutoff, maxdim)

chi = tr(cd_U_ρ_Ud)
println(imag.(chi[1])/u)
push!(char_fn,chi[1])
write_for_loop(file_name_txt, string(1), "Independent boson: T = $T, alpha = $α, N_chain = $N_chain, maxdim = $maxdim, cutoff = $cut, tau = $tau")

for t in 1:nt
    global U_ρ_Ud = apply(evol, U_ρ_Ud; cutoff, maxdim, apply_dag=true)
    global cd_U_ρ_Ud = apply(count_gates, U_ρ_Ud; cutoff, maxdim)

    chi = tr(cd_U_ρ_Ud)
    println(imag.(chi[1])/u)
    write_for_loop(file_name_txt, string(t+1), string(imag.(chi[1])/u))
    push!(char_fn,chi[1])
end
chi_pu = char_fn

mean_Q = imag.(chi_pu)/u
@show real(mean_Q)
write_to_file(file_name_txt, "Independent boson: T = $T, alpha = $α, N_chain = $N_chain, maxdim = $maxdim, cutoff = $cut, tau = $tau", string(real(mean_Q)))

p = plot()
plot!((collect(0:1:length(mean_Q)-1)*tau),real(mean_Q),label= "α = $α")
#end


plot!(legend=:bottomright)
xlabel!("t")  
title = string(L"\langle Q \rangle, N_{ch} = %$(N_chain), \beta = %$β, \mathrm{cutoff} = 10^{%$(cut)}, χ = %$maxdim")
title!(title)
display("image/png", p)
#= 

a = file_name_png;
safesave(a,p)
     =#