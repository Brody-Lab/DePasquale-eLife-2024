# DePasquale *et al* eLife 2024

Data and analysis for [DePasquale *et. al* eLife 2024](https://doi.org/10.7554/eLife.84955). 

After you clone the repository to your local machine, further steps are needed to run these notebooks. Several Julia libraries are necessary to run the analyses in each notebook. To instantiate a julia project enviorment for these notebooks and scripts, execute the following

```julia --project=. -e 'using Pkg; Pkg.instantiate()'```

This should download and add all the dependent libraries listed in `Project.toml`.

After the project has been instantiated, be sure to activate it before running any notebooks. Sometimes julia knows the current directoy is a enviorment and will activate it automatically. Otherwise, navigate to the parent directory and run `using Pkg; Pkg.activate(".")`.

Note: [PulseInputDDM.jl](https://github.com/Brody-Lab/PulseInputDDM.jl) is a custom codebase for fitting latent variable models and is required for all analyses.

## Two important notes:

- Be sure that you are using verion 0.4.4 of `PulseInputDDM`, as versions after that were significantly updated and this code will be out of date with it.
- Julia's random number generator changed after Julia 1.5, and so a small number of the notebooks here were qualitivatley (but not quantitatively) impacted by that changed.
