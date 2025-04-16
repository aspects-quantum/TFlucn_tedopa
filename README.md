## Heat fluctuations using TEDOPA 

This is a package[^SBmodel] for calculating FLUCTUATIONS of heat transfer in the Spin-Boson model[^PRX2020] using the **Time Evolving Density matrices using Orthogonal Polynomial Algorithm (_TEDOPA_)**[^Prior2010][^Chin2010]. The exact method is detailed in the paper "Heat operator approach to quantum stochastic thermodynamics in the strong-coupling regime"(https://arxiv.org/abs/2504.10631). 

We employ the **Thermofield-based chain-mapping approach for open quantum systems**[^PRA2015] that enables us to use a vacuum initial matrix product state (pure) for the environment instead of a thermal state (mixed), thereby speeding up the computation greatly.

In this package, we use the ITensor library[^Itensor] in Julia for tensor network manipulations.
&nbsp;
&NewLine;

&nbsp;

***PLEASE NOTE THE FOLLOWING POINTS before running these codes on Julia:***

We consider the following spin-boson Hamiltonian:

<img src="https://math.vercel.app/?from=%5Chat%7BH%7D%20%3D%20%5Cepsilon_0%20%5Chat%7BS%7D_z%20%2B%20%5CDelta%20%5Chat%7BS%7D_x%20%2B%20%5Csum_%7B%5Cnu%7D%20%5Comega_%7B%5Cnu%7D%20%5Chat%7Ba%7D%5E%5Cdagger_%7B%5Cnu%7D%20%5Chat%7Ba%7D_%7B%5Cnu%7D%20%2B%20%5Chat%7BS%7D_x%20%5Cotimes%20%5Csum_%7B%5Cnu%7D%20g_%7B%5Cnu%7D%20(%20%5Chat%7Ba%7D_%7B%5Cnu%7D%20%2B%20%5Chat%7Ba%7D%5E%5Cdagger_%7B%5Cnu%7D%20)" alt="Hamiltonian Equation" />

The bosonic bath has the following Ohmic spectral density:

<img src="https://math.vercel.app/?from=J(%5Comega)%3D%202%5Calpha%5Comega%5C%2C%20%5Cexp(-%5Comega%2F%5Comega_C)" alt="Spectral Density Function" />

The user can customize the above model Hamiltonian and spectral function for calculating heat fluctuations.

> Each file in **scripts** folder define the particular parameters, Hamiltonians and unitary evolutions of the considered model.

> **scripts\e0_d1_a1o25_T1_plus.jl**: independent boson model, i.e., ε = 0, Δ = 1. We take temperature, T = 1, coupling strength, α = 1.25, and initial state: |+〉.

> **scripts\e1_d0_a1o25_T0_plus.jl**: unbiased boson model, i.e., ε = 1, Δ = 0. We take temperature, T = 0, coupling strength, α = 1.25, and initial state: |+〉.

> **scripts\e1_d0_a1o25_T1_plus.jl**: unbiased boson model, i.e., ε = 1, Δ = 0. We take temperature, T = 1, coupling strength, α = 1.25, and initial state: |+〉.

&nbsp;
&NewLine;

**************************************************************************
[^SBmodel]: This code base is using the [Julia Language](https://julialang.org/) and 
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/) to make 
a reproducible scientific project named **flucn_tedopa**


[^PRX2020]: Popovic, Maria, Mark T. Mitchison, Aidan Strathearn, Brendon W. Lovett, John Goold, and Paul R. Eastham. "Quantum heat statistics with time-evolving matrix product operators." PRX Quantum 2, no. 2 (2021): 020338.

[^Prior2010]: Prior, Javier, Alex W. Chin, Susana F. Huelga, and Martin B. Plenio. "Efficient simulation of strong system-environment interactions." Physical review letters 105, no. 5 (2010): 050404.

[^Chin2010]: Chin, Alex W., Ángel Rivas, Susana F. Huelga, and Martin B. Plenio. "Exact mapping between system-reservoir quantum models and semi-infinite discrete chains using orthogonal polynomials." Journal of Mathematical Physics 51, no. 9 (2010): 092109.

[^PRA2015]: de Vega, Inés, and Mari-Carmen Banuls. "Thermofield-based chain-mapping approach for open quantum systems." Physical Review A 92, no. 5 (2015): 052116.

[^RMP2009]: Esposito, Massimiliano, Upendra Harbola, and Shaul Mukamel. "Nonequilibrium fluctuations, fluctuation theorems, and counting statistics in quantum systems." Reviews of modern physics 81, no. 4 (2009): 1665.

[^Itensor]: https://itensor.github.io/ITensors.jl/dev/index.html

## To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.


You may notice that some scripts start with the commands:
```julia
using DrWatson
@quickactivate :flucn_tedopa
```
which auto-activate the project "flucn_tedopa" and enable local path handling from DrWatson.
*****************************************************************************
