using DrWatson
@quickactivate :flucn_tedopa
## https://juliadynamics.github.io/DrWatson.jl/stable/real_world/#Making-your-project-a-usable-module-1

#@disable_warn_order
ITensors.disable_warn_order()
pgfplotsx()  # Use PGFPlotsX backend for publication-quality plots


#################            ########################################################################
################# Parameters ########################################################################
#################            ########################################################################

@show file_name_txt = string(split(split(@__FILE__, ".")[end-1], string(\))[end],".txt")

ρ0 = [1 -1;-1 1]    
ρ0 = ρ0/tr(ρ0)             ## initial spin state

ITensors.op(::OpName"ρ", ::SiteType"S=1/2") = ρ0
ITensors.op(::OpName"0", ::SiteType"Qudit", d::Int) = 1.0I[1:d,1]*1.0I[1:d,1]'
ITensors.op(::OpName"Idd", ::SiteType"Qudit", d::Int) = (1 / d) * Matrix(1.0I, d, d)
ITensors.op(::OpName"Idd", ::SiteType"S=1/2") = (1 / 2) * Matrix(1.0I, 2, 2)



cut = -11
cutoff = 10.0^cut
#maxdim_ops = 5
maxdim = 50
tau = .006             ## time step duration
nt = 600
ttotal = nt * tau           ## TOTAL TIME evolution

N_chain = 100           ## Number of chain sites for single chain-transformed environment
tot_chain = 2*N_chain+1
S_pos = N_chain+1

n1_bsn_dim = 13;
b_dim = [n1_bsn_dim-round(Int64,(n1_bsn_dim-3.6)*(i-1)/(N_chain-1))  for i = 1:N_chain]   #       ## Dimension of chain sites
boson_dim = append!(reverse(b_dim),[0],b_dim)
s_total = [(n == S_pos) ? Index(2, "S=1/2") : Index(boson_dim[n], "Qudit") for n = 1:tot_chain]
state = [(n == S_pos) ? "ρ" : "0" for n = 1:tot_chain]
ρ = normalize(MPO(s_total,state))


ω_C = 5                 ## bath cutoff
ω_0 = 0                 ## spin splitting
Ω = 1                   ## independent model if Ω = 0
if ω_0 == 1 
    model = "independent"
else
    model = "unbiased"
end

T = 5;                  ## temperature of bath
β = 1/T 
###

u = 0.001 # counting field parameter

t_list = collect(0:1:nt)*tau
#α_list = [0.1,1.5]
α = .1

support_cutoff = 10^3
supp = (0, support_cutoff)                            ## support of the weight function
Nquad = 10^7                                          ## Number of quadrature points
N_coeff = N_chain + 1
ab1 = Matrix{Float64}(undef,N_coeff,2)
ab2 = Matrix{Float64}(undef,N_coeff,2)

n(ω) = 1/(exp(β*ω) - 1)

#for α = α_list

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
    ab1[93:N_coeff, 1:2] .= repeat([a1_100 b1_100], N_coeff - 92, 1)
    ab2[93:N_coeff, 1:2] .= repeat([a2_100 b2_100], N_coeff - 92, 1)
else
    ab1 = recur_coeff(w_fn1, supp, N_coeff, Nquad)    ## recurrence coefficients for real space bath
    ab2 = recur_coeff(w_fn2, supp, N_coeff, Nquad) 
end

state = [(n == S_pos) ? "ρ" : "0" for n = 1:tot_chain]
ρ = normalize(MPO(s_total,state))

count_gates = exp_xHB(ab1, ab2, im*u, s_total)
evol = unit_gates(ω_0, Ω, c_01, c_02, ab1, ab2, tau, s_total)

char_fn=Vector{ComplexF64}()

U_ρ_Ud = ρ
cd_U_ρ_Ud = apply(count_gates, U_ρ_Ud; cutoff, maxdim)

chi = tr(cd_U_ρ_Ud)
println(imag.(chi[1])/u)
push!(char_fn,chi[1])
write_for_loop(file_name_txt, string(1), "Unbiased boson: T = $T, alpha = $α, N_chain = $N_chain, maxdim = $maxdim, cutoff = $cut, tau = $tau")

for t in 1:nt
    global U_ρ_Ud = apply(evol, U_ρ_Ud; cutoff, maxdim, apply_dag=true)
    global U_ρ_Ud = normalize(U_ρ_Ud)
    global cd_U_ρ_Ud = apply(count_gates, U_ρ_Ud; cutoff, maxdim)
    @show maxlinkdim(U_ρ_Ud)
    @show t*tau
    chi = tr(cd_U_ρ_Ud)
    println(imag.(chi[1])/u)
    write_for_loop(file_name_txt, string(t+1), string(imag.(chi[1])/u))
    push!(char_fn,chi[1])
end
chi_pu = char_fn

mean_Q = imag.(chi_pu)/u
@show real(mean_Q)
write_to_file(file_name_txt, "Unbiased boson: T = $T, alpha = $α, N_chain = $N_chain, maxdim = $maxdim, cutoff = $cut, tau = $tau", string(real(mean_Q)))

pgfplotsx()  # Use PGFPlotsX backend for publication-quality plots


default(
    fontfamily = "Times",    # Use Times New Roman for PRL
    legendfontsize = 20,     # Small legend font for clarity
    guidefontsize = 30,      # Axis labels font size
    tickfontsize = 30,       # Tick labels font size
    titlefontsize = 30,       # Title font size
    grid = true,            # Disable grid lines for a cleaner look
    linewidth = 4,         # Slightly thicker lines for visibility
    framestyle = :box,       # Box around the plot for a classic style
    size = (700, 600),       # Ensure a reasonable size for your plot
    dpi = 300,               # High resolution (300 DPI is standard for print)
    minorgrid = false,       # No minor grids
    legend = :topright       # Place the legend at the top right corner
)

p=plot((collect(0:1:length(mean_Q)-1)*tau),real(mean_Q),label= "α = $α")

p=plot((collect(0:1:length(a2)-1)*tau),real(a2),label= "|\psi_0\rangle = |-\rangle")
plot!(legend=:bottomright)
xlabel!("t")  
ylabel!(L"\langle Q \rangle")  
#= title = string(L"N_{ch} = %$(N_chain), \beta = %$β, \mathrm{cutoff} = 10^{%$(cut)}, \Delta t = %$tau, |\psi_0\rangle = |-\rangle")

 =#title!(title)
display("image/png", p)
file_name_png = string(split(split(@__FILE__, ".")[end-1], string(\))[end],"_T=$T--α=$α--$(model)--H.png")

a = file_name_png;
safesave(a,p)
#= plot(x, y,                                                                                                                                                                                                              
xlabel = "X Axis (units)",                                                                                                                                                                                          
ylabel = "Y Axis (units)",                                                                                                                                                                                          
title = "Publication-Quality Plot Example",                                                                                                                                                                         
label = L"sin(x)",                                                                                                                                                                                                  
linecolor = :blue,                                                                                                                                                                                                  
linestyle = :solid,                                                                                                                                                                                                 
marker = :circle,                                                                                                                                                                                                   
markersize = 10                                                                                                                                                                                                     
) =#
#= 
p=plot((collect(0:1:length(a1)-1)*tau),real(a1),label= L"|\psi_0\rangle = |+\rangle")
title = string(L"N_{ch} = %$(N_chain), \beta = %$β, \mathrm{cutoff} = 10^{%$(cut)}, \Delta t = %$tau")
plot!(legend=:bottomright)

xlabel!("t")  
ylabel!(L"\langle Q \rangle")  
display("image/png", p)
file_name_png = string(split(split(@__FILE__, ".")[end-1], string(\))[end],"_T=$T--α=$α--$(model).png")

a = file_name_png;
safesave(a,p) =#


