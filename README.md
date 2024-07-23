# DePasquale *et al* eLife 2023

Data and analysis for DePasquale *et. al* eLife 2023 submission.

Several Julia libraries are necessary to run the analyses in each notebook. To instantiate a julia project enviorment for these notebooks and scripts, execute the following

```julia --project=. -e 'using Pkg; Pkg.instantiate()'```

This should download and add all the dependent libraries listed in `Project.toml`.

[PulseInputDDM](https://github.com/Brody-Lab/PulseInputDDM) is a custom codebase for fitting latent variable models and is required for all analyses.