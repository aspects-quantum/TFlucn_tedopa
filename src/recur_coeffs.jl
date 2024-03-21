 #= using PolyChaos
 using Plots
 =#

function recur_coeff(w_fn, supp,N_coeff,Nquad)
    my_meas = Measure("my_meas", w_fn, supp, false, Dict())

    my_op = OrthoPoly("my_op", N_coeff-1, my_meas; Nquad)

    return coeffs(my_op)
end


#w_fn(t) = t * exp(-t);  ## weight function leading to associated Laguerre polynomials
#supp = (0, 1000)   ## support of the weight function
#Nquad = 1000;         ## Number of quadrature points
#ab = recur_coeff(w_fn,supp,N_coeff,Nquad)

