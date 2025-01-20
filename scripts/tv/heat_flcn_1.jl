using DrWatson
@quickactivate :flucn_tedopa
## https://juliadynamics.github.io/DrWatson.jl/stable/real_world/#Making-your-project-a-usable-module-1

#@disable_warn_order
ITensors.disable_warn_order()

## this script uses the function "spin_boson_time_evol" from "\src\TEBD_1.jl" ######################


#################            ########################################################################
################# Parameters ########################################################################
#################            ########################################################################

#Initial state of spin::::::::::::::::::::::
spin_state = 1 # spin "Dn"
#= spin_state = 2 # spin "Up"
spin_state = 3 # spin "X+"
spin_state = 4 # spin "X-"
=#
#::::::::::::::::::::::::::::::::::::::::::::


## time evolution parameters
cutoff = 1E-7
tau = 10^-2             ## time step duration
ttotal = 2            ## TOTAL TIME evolution

N_chain = 7            ## Number of chain sites
n1_bsn_dim = 4;
boson_dim = [n1_bsn_dim-round(Int64,(n1_bsn_dim-2)*(i-1)/(N_chain-1)) for i = 1:N_chain]           ## Dimension of chain sites


ω_C = 5                 ## bath cutoff
α = 1                   ## coupling strength
ω_0 = 1                 ## spin splitting
Ω = 0                   ## independent model if Ω = 0
c_0 = ω_C * sqrt(2 * α / pi)

T = 5;                  ## temperature of bath
β = 1/T 
###

#= N_coeff = N_chain + 1
Nquad = 5000 * N_coeff            ## Number of quadrature points
w_fn(t) = (t) * exp(-abs(t))/(1-exp(-β*t))   ## weight function leading to associated Laguerre polynomials
supp = (-10000, 10000)        ## support of the weight function

ab = recur_coeff(w_fn, supp, N_coeff, Nquad)    ## recurrence coefficients =#

ab = zeros(N_chain+1,2)
for i = (0:N_chain)
    for j = (1:2)
        if j == 1
        ab[i+1,j] = 2 * (i+1)
        else
        ab[i+1,j] = (i)*(i+1)
        end
    end
end


u = .001 # counting field 

# Make an array of 'site' INDICES for the (spin+chain)
s_total = [(n == 1) ? Index(2, "S=1/2,n=$n") : Index(boson_dim[n-1], "Qudit,n=$n") for n = 1:(N_chain+1)]
## 

chi_pu = spin_boson_char_fn(spin_state,N_chain, ω_0, Ω, ω_C, c_0, ab, s_total, tau, u, cutoff, ttotal)
chi_mu = spin_boson_char_fn(spin_state,N_chain, ω_0, Ω, ω_C, c_0, ab, s_total, tau, -u, cutoff, ttotal)
chi_0 = ones(length(chi_mu),1)

mean_Q = (-im)*(log(chi_pu)-log(chi_mu))/(2*u)

var_Q = (-1)*(log(chi_pu)+log(chi_mu))/(u^2)

# @show real(var_Q)./real(mean_Q.^2)

t_list = collect(0.0:tau:ttotal)
p=plot(t_list,real(mean_Q),label="Re[<Q>]")
plot!(t_list,imag(mean_Q),label="Im[<Q>]")
#= plot!(t_list,real(var_Q),label="Re[<<Q>>]")
plot!(t_list,imag(var_Q),label="Im[<<Q>>]") =#
#plot!(t_list,real(var_Q)./real(mean_Q.^2),label="TUR")

plot!(legend=:topright)
xlabel!("t")
display("image/png", p)
    
#=

#######################################################################################################
## plotting ###########################################################################################
#######################################################################################################

gr()

t = 0.0:tau:ttotal

x_lims = (t[1],t[end])
y_lims = (floor(minimum(Out_list)),ceil(maximum(Out_list)))

p = scatter(t, Out_list, xlims=x_lims,ylims=y_lims, label=Observable, linestyle=:dashdot,  mc=:blue, ms=5, ma=0.5,xtickfontsize=10,xlabelfontsize=12,ytickfontsize=10,
            legendfontsize=10,titlefontsize=13)  ## PLOT Out_list
plot!(legend=:bottomright)

plot!(t, Out_list,label="")
if mpo_or_not == 1
    title = string(Observable,"(t) on site = ",Obs_site,", ω_0 = ",ω_0,", Ω = ",Ω," and ",N_chain," bosons (of dim ",boson_dim,"),\nthermal initial state at β = ",β,"\n");
else
    title = string(Observable,"(t) on site = ",Obs_site,", ω_0 = ",ω_0,", Ω = ",Ω," and ",N_chain," bosons (of dim ",boson_dim, ")\n MPS initial state\n");
end
title!(title)
xlabel!("t")
display(p)

if mpo_or_not == 1
    figname = string(Observable,"_",Obs_site,"_ω=",ω_0,"_Ω=",Ω,"_N=",N_chain,"(",boson_dim,"), β=",β);
else
    figname = string(Observable,"_",Obs_site,"_ω=",ω_0,"_Ω=",Ω,"_N=",N_chain,"(", boson_dim, ")");
end
a = "plots\\$figname.png", figname;
savefig(a[1])
=#