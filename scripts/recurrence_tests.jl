using Plots, PolyChaos, LaTeXStrings, DrWatson, QuadGK


function recur_coeff(w_fn, supp,N_coeff,Nquad)
    my_meas = Measure("my_meas", w_fn, supp, false, Dict())

    my_op = OrthoPoly("my_op", N_coeff-1, my_meas; Nquad)

    return coeffs(my_op)
end

N_chain = 200           ## Number of chain sites for single chain-transformed environment

ω_C = 5                 ## bath cutoff
ω_0 = 0                 ## spin splitting
Ω = 1                   ## independent model if Ω = 0

T = 5;                  ## temperature of bath
β = 1/T 
###
α = 1

#cat(a, b) = reshape(append!(vec(a), vec(b)), size(a)[1:end-1]..., :)
support_cutoff = 200
supp = (0, support_cutoff)                            ## support of the weight function
Nquad = 10^7  
nquad = Int(log(10,Nquad))                                       ## Number of quadrature points
N_coeff = N_chain + 1
ab1 = Matrix{Float64}(undef,N_coeff,2)
ab2 = Matrix{Float64}(undef,N_coeff,2)

n(ω) = 1/(exp(β*ω) - 1)

#for α = α_list

w_fn1(k) = (2*α*k*exp(-k/ω_C))#*(1 + n(k))            ## weight function for real space bath
w_fn2(k) = (2*α*k*exp(-k/ω_C))*(n(k))                ## weight function for tilde space bath
η01 = quadgk(w_fn1, 0, support_cutoff)
c_01 = sqrt(η01[1])
η02 = quadgk(w_fn2, 0, support_cutoff)
c_02 = sqrt(η02[1])
Nmax = 92

if N_chain >= Nmax
    ab1[1:Nmax,1:2] = recur_coeff(w_fn1, supp, Nmax, Nquad)    ## recurrence coefficients for for real space bath
    ab2[1:Nmax,1:2] = recur_coeff(w_fn2, supp, Nmax, Nquad)    ## recurrence coefficients for for tilde space bath
        a1_100 = ab1[Nmax,1]; b1_100 = ab1[Nmax,2]; a2_100 = ab2[Nmax,1]; b2_100 = ab2[Nmax,2]
        [ab1[i,1] = a1_100  for i = Nmax+1:N_coeff]; [ab1[i,2] = b1_100  for i = Nmax+1:N_coeff]; 
        [ab2[i,1] = a2_100  for i = Nmax+1:N_coeff]; [ab2[i,2] = b2_100  for i = Nmax+1:N_coeff]; 
else
    ab1 = recur_coeff(w_fn1, supp, N_coeff, Nquad)    ## recurrence coefficients for real space bath
    ab2 = recur_coeff(w_fn2, supp, N_coeff, Nquad) 
end
@show sqrt(ab1[end,2])*2-ab1[end,1]
p = plot()

#= plot!(ab1[end-2:end,1],label= L"\omega_n,1")
plot!(ab2[end-2:end,1],label= L"\omega_n,2")
plot!(sqrt.(ab1[end-2:end,2]),label= L"t_n,1")
plot!(sqrt.(ab2[end-2:end,2]),label= L"t_n,2")
 =#

 plot!(ab1[:,1],label= L"\omega_n,1")
plot!(ab2[:,1],label= L"\omega_n,2")
plot!(sqrt.(ab1[:,2]),label= L"t_n,1")
plot!(sqrt.(ab2[:,2]),label= L"t_n,2")



plot!(legend=:bottomright)
xlabel!("t")  
title = string(L"N_{max} = %$(Nmax), k_{max} = %$support_cutoff, Nquad = %$nquad, \beta = %$β")
title!(title)
display("image/png", p)
#= file_name_png = string(split(split(@__FILE__, ".")[end-1], string(\))[end],".png")

a = file_name_png;
safesave(a,p)     =#
