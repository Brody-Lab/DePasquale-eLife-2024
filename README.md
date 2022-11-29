# DePasquale *et al* eLife 2022

Data and analysis for DePasquale *et. al* eLife 2022 submission. Several Julia libraries are necessary to run the analyses in each notebook. The necessary libraries are documented in each notebook. 

A Julia library can be added by entering the Julia package manager by typing `]`.

```julia
(v1.5) pkg > add [package-name]
```

[PulseInputDDM](https://github.com/Brody-Lab/PulseInputDDM) is a custom codebase for fitting latent variable models and is required for all analyses. PulseInputDDM can be added the same way:

```julia
(v1.5) pkg > add https://github.com/Brody-Lab/PulseInputDDM/
```
