using DrWatson, ITensors, ITensorMPS, PolyChaos, Plots, LinearAlgebra, QuadGK


#functions
function recur_coeff(w_fn, supp,N_coeff,Nquad)
    my_meas = Measure("my_meas", w_fn, supp, false, Dict())

    my_op = OrthoPoly("my_op", N_coeff-1, my_meas; Nquad)

    return coeffs(my_op)
end

function exp_xHB_real(ab1, x, s_tot_comp)

    #s_tot_comp = [(n == 1) ? Index(2, "S=1/2") : Index(boson_dim[n], "Qudit") for n = 1:S_pos_t]

    tot_chain = length(s_tot_comp)

    gates = ITensor[]
    i = 1
    for j in 2:tot_chain-1
        s1 = s_tot_comp[j]
        s2 = s_tot_comp[j+1]

        ω_n = ab1[i, 1]
        t_n = sqrt(ab1[i+1, 2])
        hj = ω_n * op("N", s1) * op("Id", s2) +
            t_n * op("Adag", s1) * op("A", s2) +
            t_n * op("A", s1) * op("Adag", s2)
        
        Gj = exp((x / 2)* hj)
        i=i+1
        push!(gates, Gj)
    end
    append!(gates, reverse(gates))

    return gates
end

##########################################################################
function HBB_tf(ab1, ab2, s_total)

    #s_total = [(n == S_pos_t) | (n == S_pos_r) ? Index(2, "S=1/2") : Index(boson_dim, "Qudit") for n = 1:tot_chain]
    tot_chain = length(s_total)
    S_pos_t = Int(tot_chain/2)
    S_pos_r = S_pos_t+1

    hbb = OpSum()
    i=1
    for j in S_pos_r+1:tot_chain-1
        ω_n = ab1[i, 1]
        t_n = sqrt(ab1[i+1, 2])
        hbb += ω_n, "N", j
        hbb += t_n, "Adag", j, "A", j+1
        hbb += t_n, "A", j, "Adag", j+1
        i = i+1
    end

    i=1
    for j in S_pos_t-1:-1:2
        s1 = s_total[j]
        s2 = s_total[j-1]
        ω_n = ab2[i, 1]
        t_n = sqrt(ab2[i+1, 2])
        hbb += -ω_n, "N", j
        hbb += -t_n, "Adag", j, "A", j-1
        hbb += -t_n, "A", j, "Adag", j-1
        i = i+1
    end
    
    HBB = MPO(hbb, s_total)

    return HBB
end


ITensors.disable_warn_order()

#################            ########################################################################
################# Parameters ########################################################################
#################            ########################################################################


ρ = [1 0;0 0]    
ρ = ρ/tr(ρ)             ## initial spin state

ITensors.op(::OpName"ρ", ::SiteType"S=1/2") = ρ


cutoff = 1E-9
#maxdim_ops = 5
maxdim = 10
tau = .01             ## time step duration
nt = 6
Nbeta = 10
ttotal = nt * tau           ## TOTAL TIME evolution

N_chain = 15            ## Number of chain sites for single chain-transformed environment
tot_chain = 2*N_chain+2


# Make an array of 'site' INDICES for the (spin+chain)
S_pos_t = N_chain+1
S_pos_r = N_chain+2

println(S_pos_r)

n1_bsn_dim = 7;
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
α = .5

#cat(a, b) = reshape(append!(vec(a), vec(b)), size(a)[1:end-1]..., :)
support_cutoff = 7*10^2
supp = (0, support_cutoff)                            ## support of the weight function
Nquad = 10^7                                          ## Number of quadrature points
N_coeff = N_chain + 1

ab1 = Matrix{Float64}(undef,N_coeff,2)
ab2 = Matrix{Float64}(undef,N_coeff,2)

n(k) = 1/(exp(β*k) - 1)

w_fn1(k) = (2*α*k*exp(-k/ω_C))*(1 + n(k))           ## weight function for real space bath
w_fn2(k) = (2*α*k*exp(-k/ω_C))*n(k)           ## weight function for tilde space bath
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
    ab1 = recur_coeff(w_fn1, supp, N_coeff, Nquad)    ## recurrence coefficients for for real space bath
    ab2 = recur_coeff(w_fn2, supp, N_coeff, Nquad)    ## recurrence coefficients for for tilde space bath
end


maxdim_init = 20

tot_chain = length(s_total)

hbb = HBB_tf(ab1, ab2, s_total)
hbb_comp = MPO(s_tot_comp)
for i = 1:Int(length(s_tot_comp))
    hbb_comp[S_pos_t+1-i] = hbb[i]*hbb[tot_chain+1-i]
end



#= 
r0 = apply(exp_xHB_tf(ab1, -0.5*β/Nbeta, s_total), MPO(s_total, ["Idd" for i=1:tot_chain]); cutoff, maxdim = maxdim_init)
r0_comp = MPO(s_tot_comp)
for i = 1:Int(length(s_tot_comp))
    r0_comp[S_pos_t+1-i] = r0[i]*r0[tot_chain+1-i]
end =#

#r0 = apply(exp_xHB_real(ab1, -0.5*β/Nbeta, s_tot_comp), MPO(s_tot_comp, ["Idd" for i=1:length(s_tot_comp)]); cutoff, maxdim = maxdim_init)
I0 = MPO(s_tot_comp, "Id")
I0_mps = convert(MPS, I0)

for i = 1:Int(length(s_tot_comp))
    I0_mps[i] = I0_mps[i]*delta(s_tot_comp[i]', s_total[S_pos_t+1-i])
end

r0_mps = tdvp(hbb_comp, -u*1.0im, I0_mps; nsteps=5, updater_backend="exponentiate", cutoff=1e-9, maxdim, normalize=true)

lv = left_vac(s_total)
S_pos_t = Int(tot_chain/2)
S_pos_r = S_pos_t + 1
for i = 1:S_pos_t
    r0[i] = lv[i]*r0[i]
end
rh0 = MPS(s_total)
for i = 1:S_pos_t
    if i == 1
        rh0[1], rh0[tot_chain] = qr(r0[1]*r0[tot_chain], (s_total[1]', commoninds(r0[1],r0[2])[1]))
    elseif i == S_pos_t
        rh0[S_pos_t], rh0[S_pos_r] = qr(r0[S_pos_t]*r0[S_pos_r], (s_total[S_pos_t]', commoninds(r0[S_pos_t],r0[S_pos_t-1])[1]))
    else
        rh0[Int(i)], rh0[Int(tot_chain-i+1)] = qr(r0[Int(i)]*r0[Int(tot_chain-i+1)], (s_total[Int(i)]', commoninds(r0[Int(i)],r0[Int(i-1)])[1], commoninds(r0[Int(i)],r0[Int(i+1)])[1]))
    end
end
normalize!(rh0)
rh0 = replaceprime(rh0, 1=>0)

for i = 1:Nbeta-1           ## perform IMAGINARY TIME evolution
    rh0 = apply(exp_xHB_tf(ab1, -0.5*β/Nbeta, s_total), rh0; cutoff, maxdim = maxdim)
    normalize!(rh0)
    println("thermal state", i)
end