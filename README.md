# Gauge Monte Carlo

This is a proof-of-concept implementation of a 1D canonical Monte Carlo simulation with phase space shifting as described in:

**Why gauge invariance applies to statistical mechanics**  
*Johanna Müller, Florian Sammüller, and Matthias Schmidt, [J. Phys. A: Math. Theor. (to be published)](https://doi.org/10.1088/1751-8121/adbfe6) (2024); [arXiv:2409.14166](https://arxiv.org/abs/2409.14166).*


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
