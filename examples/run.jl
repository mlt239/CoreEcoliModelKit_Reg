using CoreEcoliModel

# setup calculation -
path_to_measurements_file = "./Emmerling-JBac-2002-MCM.json"

# estimate the flux -
(flux_distribution, results_array) = maximize_specific_growth_rate(path_to_measurements_file)

# compute some flux ratios (using the mean, we should really be sampling these ...)
qGlucose = flux_distribution[50,1]
scaledFluxArray = (1/qGlucose)*flux_distribution[:,1]

# Look at the PPP split -
PGI_flux = 100*scaledFluxArray[74]
PPP_flux = 100*scaledFluxArray[48]

# A-fluxes
PEP_OAA_flux = 100*scaledFluxArray[79] # PPC
OAA_PEP_flux = 100*scaledFluxArray[80] # PPCK
ME1_flux = 100*scaledFluxArray[65]  # ME1 (nadh)
ME2_flux = 100*scaledFluxArray[66]  # ME1 (nadph)
