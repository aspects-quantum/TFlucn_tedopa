## Heat fluctuations using TEDOPA 

This is a package[^SBmodel] for calculating FLUCTUATIONS of heat transfer in the Spin-Boson model[^PRX2020] using the **Time Evolving Density matrices using Orthogonal Polynomial Algorithm (_TEDOPA_)**[^Prior2010][^Chin2010]. 

We employ the **Thermofield-based chain-mapping approach for open quantum systems**[^PRA2015] that enables us to use a vacuum initial matrix product state (pure) for the environment instead of a thermal state (mixed), thereby speeding up the computation greatly.

In this package, we use the ITensor library[^Itensor] in Julia for tensor network manipulations.
&nbsp;
&NewLine;

&nbsp;

***PLEASE NOTE THE FOLLOWING POINTS before running these codes on Julia:***

> The **scripts** files define all the parameters of the model: **scripts\mean_Q.jl** (computes mean heat exchange), **scripts\var_Q.jl** (computes variance of heat exchange), **scripts\obs_exp.jl** (computes expected value of any observable).

> Hamiltonians and Gates are defined in the file **src\ham_gates.jl**.

> Initial state of the (system+environment) is defined in the file **src\init_state.jl**.

> The heat characteristic function[^RMP2009] is calculated in the file **src\charfn.jl**.

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
