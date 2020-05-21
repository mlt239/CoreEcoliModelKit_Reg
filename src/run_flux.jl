#This is what you should call to actually perform a regulated and non regulated flux

function run_flux()
#using CoreEcoliModel

# change this
path_to_measurements_file = "C:\\Users\\matt-\\Box\\Fischbach Research\\FBA stuff\\CoreEcoliModelKit\\examples\\Emmerling-JBac-2002-MCM.json"

# estimate the flux -
results_array = maximize_specific_growth_rate_regulated(path_to_measurements_file)
(results_array2, fluxes) = maximize_specific_growth_rate_non_regulated(path_to_measurements_file)
# compute some flux ratios (using the mean, we should really be sampling these ...)
#=qGlucose = flux_distribution[50,1]
scaledFluxArray = (1/qGlucose)*flux_distribution[:,1]

# Look at the PPP split -
PGI_flux = 100*scaledFluxArray[74]
PPP_flux = 100*scaledFluxArray[48]

# A-fluxes
PEP_OAA_flux = 100*scaledFluxArray[79] # PPC
OAA_PEP_flux = 100*scaledFluxArray[80] # PPCK
ME1_flux = 100*scaledFluxArray[65]  # ME1 (nadh)
ME2_flux = 100*scaledFluxArray[66]  # ME1 (nadph) =#
return (results_array, results_array2, fluxes)
end

function run_flux(initialReg)
#using CoreEcoliModel

# setup calculation -
path_to_measurements_file = "C:\\Users\\matt-\\Box\\Fischbach Research\\FBA stuff\\CoreEcoliModelKit\\examples\\Emmerling-JBac-2002-MCM.json"

# estimate the flux -
results_array = maximize_specific_growth_rate_regulated(path_to_measurements_file, initialReg)

# compute some flux ratios (using the mean, we should really be sampling these ...)
#=qGlucose = flux_distribution[50,1]
scaledFluxArray = (1/qGlucose)*flux_distribution[:,1]

# Look at the PPP split -
PGI_flux = 100*scaledFluxArray[74]
PPP_flux = 100*scaledFluxArray[48]

# A-fluxes
PEP_OAA_flux = 100*scaledFluxArray[79] # PPC
OAA_PEP_flux = 100*scaledFluxArray[80] # PPCK
ME1_flux = 100*scaledFluxArray[65]  # ME1 (nadh)
ME2_flux = 100*scaledFluxArray[66]  # ME1 (nadph) =#
end
