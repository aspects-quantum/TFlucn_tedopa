using DrWatson
@quickactivate :flucn_tedopa
## https://juliadynamics.github.io/DrWatson.jl/stable/real_world/#Making-your-project-a-usable-module-1

#@disable_warn_order
ITensors.disable_warn_order()

#################            ########################################################################
################# Parameters ########################################################################
#################            ########################################################################

@show file_name_txt = string(split(split(@__FILE__, ".")[end-1], string(\))[end],".txt")

ρ0 = [1 1;1 1]    
ρ0 = ρ0/tr(ρ0)             ## initial spin state

ITensors.op(::OpName"ρ", ::SiteType"S=1/2") = ρ0



cut = -11
cutoff = 10.0^cut
#maxdim_ops = 5
maxdim = 1000
tau = .001             ## time step duration
nt = 500
ttotal = nt * tau           ## TOTAL TIME evolution

N_chain = 100           ## Number of chain sites for single chain-transformed environment
tot_chain = 2*N_chain+1
S_pos = N_chain+1

println(S_pos)

n1_bsn_dim = 7;
b_dim = [n1_bsn_dim-round(Int64,(n1_bsn_dim-2.6)*(i-1)/(N_chain-1))  for i = 1:N_chain]   #       ## Dimension of chain sites
boson_dim = append!(reverse(b_dim),[0],b_dim)

s_total = [(n == S_pos) ? Index(2, "S=1/2") : Index(boson_dim[n], "Qudit") for n = 1:tot_chain]

ITensors.op(::OpName"0", ::SiteType"Qudit", d::Int) = 1.0I[1:d,1]*1.0I[1:d,1]'
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


u = 0.001 # counting field parameter

t_list = collect(0:1:nt)*tau
#α_list = [0.1,1.5]
α = .1

#cat(a, b) = reshape(append!(vec(a), vec(b)), size(a)[1:end-1]..., :)
support_cutoff = 700
supp = (0, support_cutoff)                            ## support of the weight function
Nquad = 10^7                                          ## Number of quadrature points
N_coeff = N_chain + 1
ab1 = Matrix{Float64}(undef,N_coeff,2)
ab2 = Matrix{Float64}(undef,N_coeff,2)


T = .1;                  ## temperature of bath
β = 1/T 
###
n(ω) = 1/(exp(β*ω) - 1)

#for α = α_list
#= 
w_fn1(k) = (2*α*k*exp(-k/ω_C))*(1 + n(k))            ## weight function for real space bath
w_fn2(k) = (2*α*k*exp(-k/ω_C))*(n(k))                ## weight function for tilde space bath
η01 = quadgk(w_fn1, 0, support_cutoff)
c_01 = sqrt(η01[1])
η02 = quadgk(w_fn2, 0, support_cutoff)
c_02 = sqrt(η02[1])

if N_chain >= 92
    ab1[1:92,1:2] = recur_coeff(w_fn1, supp, 92, Nquad)    ## recurrence coefficients for for real space bath
    ab2[1:92,1:2] = recur_coeff(w_fn2, supp, 92, Nquad)    ## recurrence coefficients for for tilde space bath
        a1_100 = ab1[92,1]; b1_100 = ab1[92,2]; a2_100 = ab2[92,1]; b2_100 = ab2[92,2]
        [ab1[i,1] = a1_100  for i = 93:N_coeff]; [ab1[i,2] = b1_100  for i = 93:N_coeff]; 
        [ab2[i,1] = a2_100  for i = 93:N_coeff]; [ab2[i,2] = b2_100  for i = 93:N_coeff]; 
else
    ab1 = recur_coeff(w_fn1, supp, N_coeff, Nquad)    ## recurrence coefficients for real space bath
    ab2 = recur_coeff(w_fn2, supp, N_coeff, Nquad) 
end =#
#= 

# Saving variables
save("ab.jld2", "ab1", ab1, "ab2", ab2)

# Loading variables
loaded_vars = load("ab.jld2")
ab1 = loaded_vars["ab1"]
ab2 = loaded_vars["ab2"]
 =#
state = [(n == S_pos) ? "ρ" : "0" for n = 1:tot_chain]
ρ = normalize(MPO(s_total,state))

count_gates1 = exp_xHB(ab1, ab2, im*u, s_total)
#count_gates2 = exp_xHB(ab1, ab2, -im*u, s_total)
evol = unit_gates(ω_0, Ω, c_01, c_02, ab1, ab2, tau, s_total)

char_fn_1=Vector{ComplexF64}()
char_fn_2=Vector{ComplexF64}()


U_ρ_Ud1 = ρ
cd_U_ρ_Ud1 = apply(count_gates1, U_ρ_Ud1; cutoff, maxdim)
#= global U_ρ_Ud2 = ρ
cd_U_ρ_Ud2 = apply(count_gates2, ρ; cutoff, maxdim)
 =#
@show chi_1 = tr(cd_U_ρ_Ud1)
@show chi_2 = conj(chi_1)

#println(real(-(log.(chi_2)+log.(chi_1))/(u^2)))
println(real(-(4*(chi_1+chi_2-2)-(chi_1-chi_2)^2)/(4*u^2)))  
push!(char_fn_1,chi_1[1])
push!(char_fn_2,chi_2[1])

write_for_loop(file_name_txt, string(1), "$(model) boson: T = $T, alpha = $α, N_chain = $N_chain, maxdim = $maxdim, cutoff = $cut, tau = $tau")
for t in 1:nt
    global U_ρ_Ud1 = apply(evol, U_ρ_Ud1; cutoff, maxdim, apply_dag=true)
    global U_ρ_Ud1 = normalize(U_ρ_Ud1)
    global cd_U_ρ_Ud1 = apply(count_gates1, U_ρ_Ud1; cutoff, maxdim)
     chi_1 = tr(cd_U_ρ_Ud1)
     chi_2 = conj(chi_1)
    #println(real(-(log.(chi_2[1])-2*log.(chi_1[1]))/(u^2)))
    #println(real(-(log.(chi_2[1])+log.(chi_1[1]))/(u^2)))
    println(real(-(4*(chi_1+chi_2-2)-(chi_1-chi_2)^2)/(4*u^2)))
    push!(char_fn_1,chi_1[1])
    push!(char_fn_2,chi_2[1])
    write_for_loop(file_name_txt, string(t+1), string(real(-(log.(chi_2[1])+log.(chi_1[1]))/(u^2))))
end
#var_Q = real(-(log.(char_fn_2)+log.(char_fn_1))/(u^2))
var_Q = real(-(4.0.*(char_fn_1 + char_fn_2 - 2.0.*ones(length(char_fn_1),1)) - (char_fn_1 - char_fn_2).^2)/(4*u^2))
write_to_file(file_name_txt, "$(model) boson: T = $T, alpha = $α, N_chain = $N_chain, maxdim = $maxdim, cutoff = $cut, tau = $tau", string(var_Q))

p = plot()
plot!((collect(0:1:length(var_Q)-1)*tau)[begin:10:end],var_Q[begin:10:end],label= "α = $α")
#end

plot!(legend=:bottomright)
xlabel!("t")  
title = string(L"\langle Q^2 \rangle, N_{ch} = %$(N_chain), \beta = %$β, \mathrm{cutoff} = 10^{%$(cut)}, \Delta t = %$tau")
title!(title)
display("image/png", p)
file_name_png = string(split(split(@__FILE__, ".")[end-1], string(\))[end],"_T=$T--α=$α--$(model).png")

a = file_name_png;
safesave(a,p)
    
