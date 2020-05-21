module CoreEcoliModel

# internal includes -
include("Include.jl")

#
# maximize_specific_growth_rate(path_to_measurements_file::String; number_of_samples::Int64 = 100)
export maximize_specific_growth_rate
export maximize_flux_at_index
export moma_calculation
export flux_variabilty_analysis

# Data dictionary methods -
export generate_default_data_dictionary
export constrain_measured_fluxes
export constrain_measured_metabolites
export optimize_specific_growth_rate

# Export utility methods -
export write_model_reactions_to_file

# export types -
export VLOptimalFluxResult

end # module
