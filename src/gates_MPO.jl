using ITensors

N = 20
s = SpinHalf(N, ConserveQNs=false)

tau = 0.1
G = Vector{ITensor}(undef, N)
for j in 1:N-1
    G[j] = expHermitian(op(s, "Sx", j)*op(s, "Sx", j+1), -tau)
end


#Split up the gates
A = Vector{ITensor}(undef, N)
B = Vector{ITensor}(undef, N)
for j in 1:N-1
    Aj, Bj = factor(G[j], [s(j), prime(s(j))])
    A[j] = Aj
    B[j] = Bj
end

#Some temporary internal indices
t = Vector{Index}(undef, N+1)
for j in 1:N
    t[j] = sim(s(j))
end

#Replace internal indices so that the correct indices contract
for j in 2:N-1
    B[j-1] *= delta(prime(s(j)), t[j])
    A[j] *= delta(s(j), t[j])
end

#Create and store the MPO tensors
M = MPO(N)
M[1] = A[1]
for j in 2:N-1
    M[j] = B[j-1]*A[j]
end
M[N] = B[N-1]

#As a test, apply this MPO to a state and compare to applying the gates

init = InitState(s)
for j in 1:N
    init[j] = "Up"
end
psi = MPS(init)

Mpsi = applyMPO(M, psi)

Gpsi = copy(psi)

position(Gpsi, 1)
for j in 1:N-1
    phi = Gpsi[j]*Gpsi[j+1]*G[j]
    noPrime!(phi)
    U, S, V = svd(phi, inds(Gpsi[j]))
    Gpsi[j] = U
    Gpsi[j+1] = S*V
end

println(inner(Mpsi, Mpsi))
println(inner(Gpsi, Gpsi))
println(inner(Gpsi, Mpsi))

