# # includes -
# include("Include.jl")
#
# # setup some global consts -
# const path_to_cobra_mat_file = "$(pwd())/cobra/config/matlab_cobra_files/modelReg.mat"
# const model_file_name = "modelReg"
# const organism_id = :eco    # we use the KEGG organism symbols
# const path_to_measurements_file = "$(pwd())/experimental_data/test_data/Glucose.json"

function moma_calculation(flux_soln_array, flux_bounds_array, ko_index)

    # load the default data_dictionary -
    default_data_dictionary = generate_default_data_dictionary(path_to_cobra_mat_file, model_file_name, organism_id);

    # pass the default dictionary to a customization method -
    updated_data_dictionary = optimize_specific_growth_rate(default_data_dictionary);

    # KO -
    flux_bounds_array[ko_index,1] = 0.0
    flux_bounds_array[ko_index,2] = 0.0
    updated_data_dictionary["flux_bounds_array"] = flux_bounds_array

    # execute MOMA calculation -
    moma_soln = moma(flux_soln_array, updated_data_dictionary)

    # return -
    return moma_soln
end

# # call an executable function -
# (flux_distribution, results_array) = maximize_specific_growth_rate(path_to_measurements_file; number_of_samples = 1000);
#
# # take the mean -
# wild_type_soln = results_array[10].flux_array
# wild_type_soln_fba = results_array[10].flux_bounds_array
# moma_soln = moma_calculation(wild_type_soln, wild_type_soln_fba, 12)
