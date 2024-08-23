module DePasqualeElife2024

using Flatten, Parameters, Distributions, PulseInputDDM, Distributed
import PulseInputDDM: nθparams, optimize, diffLR, DDM, DDMθ, neuraldata
import PulseInputDDM: initialize_latent_model, choice_likelihood!, adapt_clicks, latent_one_step!
import StatsFuns: logistic
import Optim: converged, minimizer

export neural_choice_GLM_DDM, loglikelihood, θneural_choice_GLM

include("neural-choice_GLM_model.jl")
include("neural-choice_model-sep.jl")
include("neural_model-sep.jl")
include("xcorrs.jl")
include("change-of-mind.jl")

end