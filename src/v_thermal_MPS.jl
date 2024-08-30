function thermal_MPS(ab1, s_total, β, Nbeta, cutoff, maxdim)
    tot_chain = length(s_total)

    maxdim_init = 20
    r0 = apply(exp_xHB_tf(ab1, -0.5*β/Nbeta, s_total), MPO(s_total, ["Idd" for i=1:tot_chain]); cutoff, maxdim = maxdim_init)
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

    for i = 2:Nbeta           ## perform IMAGINARY TIME evolution
        rh0 = apply(exp_xHB_tf(ab1, -0.5*β/Nbeta, s_total), rh0; cutoff, maxdim = maxdim)
        normalize!(rh0)
        println("thermal state", i)
    end
    #= f = h5open("myfile.h5","w")
    write(f,"thermal_state", rh0)
    close(f) =#

    return rh0
end


function left_vac(s_total)
    tot_chain = length(s_total)
    S_pos_t = Int(tot_chain/2)
    vec = ITensor[]
    for i=1:S_pos_t
        push!(vec, delta(s_total[Int(i)],s_total[Int(tot_chain-i+1)]))
    end
    return vec
end

#= 
```julia
│ s = siteinds("S=1/2")
│ psi = randomMPS(s)
│ H = MPO(s, "Id")
│ inner(psi, H, psi)
│ ```
│ 
│ `psi` has the Index structure `-s-(psi)` and `H` has the Index structure
│ `-s'-(H)-s-`, so the Index structure of would be `(dag(psi)-s- -s'-(H)-s-(psi)`
│  unless the prime levels were fixed. Previously we tried fixing the prime level
│   in situations like this, but we will no longer be doing that going forward.
│ 
│ There are a few ways to fix this. You can simply change:
│ 
│ ```julia
│ inner(psi, H, psi)
│ ```
│ 
│ to:
│ 
│ ```julia
│ inner(psi', H, psi)
│ ```
│ 
│ in which case the Index structure will be `(dag(psi)-s'-(H)-s-(psi)`.
│ 
│ Alternatively, you can use the `Apply` function:
│ 
│ ```julia
│ 
│ inner(psi, Apply(H, psi))
│ ```
│ 
│ In this case, `Apply(H, psi)` represents the "lazy" evaluation of
│ `apply(H, psi)`. The function `apply(H, psi)` performs the contraction of
│ `H` with `psi` and then unprimes the results, so this versions ensures that
│ the prime levels of the inner product will match.
│ 
│ Although the new behavior seems less convenient, it makes it easier to
│ generalize `inner(::MPS, ::MPO, ::MPS)` to other types of inputs, like `MPS`
│ and `MPO` with different tag and prime conventions, multiple sites per tensor,
│ `ITensor` inputs, etc.
│  =#