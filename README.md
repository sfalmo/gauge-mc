# Gauge Monte Carlo

This is a proof-of-concept implementation of a 1D canonical Monte Carlo simulation with phase space shifting as described in:

**Why gauge invariance applies to statistical mechanics**  
*Johanna Müller, Florian Sammüller, and Matthias Schmidt, [J. Phys. A: Math. Theor. (to be published)](https://doi.org/10.1088/1751-8121/adbfe6) (2025); [arXiv:2409.14166](https://arxiv.org/abs/2409.14166).*

See also:

**Gauge invariance of equilibrium statistical mechanics**  
*Johanna Müller, Sophie Hermann, Florian Sammüller, and Matthias Schmidt, [Phys. Rev. Lett. **133**, 217101](https://doi.org/10.1103/PhysRevLett.133.217101) (2024); [arXiv:2406.19235](https://arxiv.org/abs/2406.19235).*


## Instructions

A recent version of [Julia](https://julialang.org/downloads/) needs to be installed on your system.
Launch the Julia interpreter within this directory and type `]` to enter the package manager.
Activate the environment and install the required packages as follows:

```julia
activate .
instantiate
```

Type backspace to exit the package manager.
Run the main program:

```julia
include("main.jl")
```

The program performs Monte Carlo simulations with and without phase space shifting and plots the results, which will take a few minutes.
