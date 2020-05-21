# # includes -
# include("Include.jl")
#
# # setup some global consts -
# const path_to_cobra_mat_file = "$(pwd())/cobra/config/matlab_cobra_files/modelReg.mat"
# const model_file_name = "modelReg"
# const organism_id = :eco    # we use the KEGG organism symbols
# const path_to_measurements_file = "$(pwd())/experimental_data/test_data/Glucose.json"
#Added new function for regulated FBA

"""
TODO: Fill me in with some stuff ...
"""
function maximize_flux_at_index(flux_index_array::Array{Int64,1}, path_to_measurements_file::String; number_of_samples::Int64 = 100)

    # initalize -
    results_array = Array{VLOptimalFluxResult,1}()

    # Declare a progress meter for user feedback -
    p = Progress(number_of_samples,color=:yellow)

    # ok, so lets sample ...
    for sample_index = 1:number_of_samples

        # load the default data_dictionary -
        default_data_dictionary = generate_default_data_dictionary();

        # pass the default dictionary to a customization method -
        updated_data_dictionary = optimize_flux_at_index(flux_index_array, default_data_dictionary);

        # update dictionary with experimental data?
        updated_data_dictionary_fluxes = constrain_measured_fluxes(updated_data_dictionary, path_to_measurements_file);
        updated_data_dictionary_metabolites = constrain_measured_metabolites(updated_data_dictionary_fluxes, path_to_measurements_file);

        # estimate the optimal flux distrubution -
        (objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag) = calculate_optimal_flux_distribution(updated_data_dictionary_metabolites);

        # build a return type -
        fluxResult = VLOptimalFluxResult();
        fluxResult.objective_value = objective_value;
        fluxResult.flux_array = calculated_flux_array;
        fluxResult.dual_array = dual_value_array;
        fluxResult.uptake_array = uptake_array;
        fluxResult.exit_flag = exit_flag;
        fluxResult.status_flag = status_flag;

        # get some problem setup information -
        fluxResult.flux_bounds_array = updated_data_dictionary["flux_bounds_array"];

        # cache -
        push!(results_array, fluxResult);

        # user message -
        msg = "Completed $(sample_index) of $(number_of_samples) trials ...";

        # update the progress bar -
        ProgressMeter.next!(p; showvalues = [(:status,msg)]);
    end

    @info "Completed ...\r";

    # how many flux are there?
    number_of_fluxes = length(results_array[1].flux_array);

    # compute the flux values -
    flux_ensemble = zeros(number_of_fluxes,1);
    for flux_object in results_array

        # grab the flux -
        flux_array = flux_object.flux_array;

        # cache -
        flux_ensemble = [flux_ensemble flux_array];
    end

    # cut off the zeros -
    flux_ensemble = flux_ensemble[:,2:end];

    # compute the mean and std -
    µ = mean(flux_ensemble, dims=2);
    σ = std(flux_ensemble, dims=2);
    flux_distribution = [µ σ];

    # return -
    return (flux_distribution, results_array)
end

"""
TODO: Fill me in with some stuff ...
"""
function maximize_specific_growth_rate(path_to_measurements_file::String; number_of_samples::Int64 = 100)

    # initalize -
    results_array = Array{VLOptimalFluxResult,1}()

    # Declare a progress meter for user feedback -
    p = Progress(number_of_samples,color=:yellow)

    # ok, so lets sample ...
    for sample_index = 1:number_of_samples

        # load the default data_dictionary -
        default_data_dictionary = generate_default_data_dictionary();

        # pass the default dictionary to a customization method -
        updated_data_dictionary = optimize_specific_growth_rate(default_data_dictionary);

        # update dictionary with experimental data?
        updated_data_dictionary = constrain_measured_fluxes(updated_data_dictionary, path_to_measurements_file);

        # estimate the optimal flux distrubution -
        (objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag) = calculate_optimal_flux_distribution(updated_data_dictionary);

        # build a return type -
        fluxResult = VLOptimalFluxResult();
        fluxResult.objective_value = objective_value;
        fluxResult.flux_array = calculated_flux_array;
        fluxResult.dual_array = dual_value_array;
        fluxResult.uptake_array = uptake_array;
        fluxResult.exit_flag = exit_flag;
        fluxResult.status_flag = status_flag;

        # get some problem setup information -
        fluxResult.flux_bounds_array = updated_data_dictionary["flux_bounds_array"];

        # cache -
        push!(results_array, fluxResult);

        # user message -
        msg = "Completed $(sample_index) of $(number_of_samples) trials ...";

        # update the progress bar -
        ProgressMeter.next!(p; showvalues = [(:status,msg)]);
    end

    @info "Completed ...\r";

    # how many flux are there?
    number_of_fluxes = length(results_array[1].flux_array);

    # compute the flux values -
    flux_ensemble = zeros(number_of_fluxes,1);
    for flux_object in results_array

        # grab the flux -
        flux_array = flux_object.flux_array;

        # cache -
        flux_ensemble = [flux_ensemble flux_array];
    end

    # cut off the zeros -
    flux_ensemble = flux_ensemble[:,2:end];

    # compute the mean and std -
    µ = mean(flux_ensemble, dims=2);
    σ = (1/sqrt(number_of_samples))*std(flux_ensemble, dims=2);
    flux_distribution = [µ σ];

    # return -
    return (flux_distribution, results_array)
end

function sample_flux_space_with_experimental_constraints(solution_bounds_array::Array{Float64,2}, number_of_samples::Int64)

    # load the default data_dictionary -
    default_data_dictionary = generate_default_data_dictionary(path_to_cobra_mat_file, model_file_name, organism_id);

    # update dictionary with experimental data?
    updated_data_dictionary = constrain_measured_fluxes(default_data_dictionary, path_to_measurements_file)

    # call the sample method -
    return sample_flux_space(solution_bounds_array, updated_data_dictionary,number_of_samples)
end

function maximize_specific_growth_rate_regulated(path_to_measurements_file::String, initialReg)

    # initalize -
    results_array = Array{VLOptimalFluxResult,1}()



    # ok, so lets sample ...


        # load the default data_dictionary -
        default_data_dictionary = generate_default_data_dictionary();

        # pass the default dictionary to a customization method -
        updated_data_dictionary = optimize_specific_growth_rate(default_data_dictionary);

        # update dictionary with experimental data?
        updated_data_dictionary = constrain_measured_fluxes(updated_data_dictionary, path_to_measurements_file);


        # estimate the optimal flux distrubution -
        #(objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag) = calculate_optimal_flux_distribution(updated_data_dictionary);
        (regulated_data_dictionary, calculated_flux_array, dual_value_array, uptake_array, DRgenes, constrainedRxns, states, exit_flag, status_flag) = optimizeRegModel(updated_data_dictionary, initialReg)

        # build a return type -
        #fluxResult = VLOptimalFluxResult();
        fluxResult.objective_value = objective_value;
        fluxResult.flux_array = calculated_flux_array;
        fluxResult.dual_array = dual_value_array;
        fluxResult.uptake_array = uptake_array;
        fluxResult.exit_flag = exit_flag;
        fluxResult.status_flag = status_flag;
        fluxResults.DRgenes = DRgenes;
        fluxResults.constrainedRxns = constrainedRxns;
        fluxResults.states = states;

        # get some problem setup information -
        fluxResult.flux_bounds_array = regulated_data_dictionary["flux_bounds_array"];

        # cache -
        push!(results_array, fluxResult);

    # how many flux are there?
    #number_of_fluxes = length(results_array[1].flux_array);

    # compute the flux values -
    #flux_ensemble = zeros(number_of_fluxes,1);
    #for flux_object in results_array

        # grab the flux -
        #flux_array = flux_object.flux_array;

        # cache -
        #flux_ensemble = [flux_ensemble flux_array];
    #end

    # cut off the zeros -
    #flux_ensemble = flux_ensemble[:,2:end];

    # compute the mean and std -
    #µ = mean(flux_ensemble, dims=2);
    #σ = (1/sqrt(number_of_samples))*std(flux_ensemble, dims=2);
    #flux_distribution = [µ σ];

    # return -
    return (results_array)
end

#The 2 functions below are shortcuts to performing regulated and non-regulated FBA
function maximize_specific_growth_rate_regulated(path_to_measurements_file::String;) #using optimizeRegModel

    # initalize -
    results_array = Array{VLOptimalFluxResult,1}()



    # ok, so lets sample ...


        # load the default data_dictionary -
        default_data_dictionary = generate_default_data_dictionary();

        # pass the default dictionary to a customization method -
        updated_data_dictionary = optimize_specific_growth_rate(default_data_dictionary);

        # update dictionary with experimental data?
        updated_data_dictionary = constrain_measured_fluxes(updated_data_dictionary, path_to_measurements_file);
        #updated_data_dictionary["flux_bounds_array"][83,:] .= 0.0

        # estimate the optimal flux distrubution -
        #(objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag) = calculate_optimal_flux_distribution(updated_data_dictionary);
        results_array = optimizeRegModel(updated_data_dictionary)

        # build a return type -
        #fluxResult = VLOptimalFluxResult();
        #fluxResult.objective_value = objective_value;
        #fluxResult.flux_array = calculated_flux_array;
        #fluxResult.dual_array = dual_value_array;
        #fluxResult.uptake_array = uptake_array;
        #fluxResult.exit_flag = exit_flag;
        #fluxResult.status_flag = status_flag;
        #fluxResults.DRgenes = DRgenes;
        #fluxResults.constrainedRxns = constrainedRxns;
        #fluxResults.states = states;

        # get some problem setup information -
        #fluxResult.flux_bounds_array = regulated_data_dictionary["flux_bounds_array"];

        # cache -
        #push!(results_array, fluxResult);

    # how many flux are there?
    #number_of_fluxes = length(results_array[1].flux_array);

    # compute the flux values -
    #flux_ensemble = zeros(number_of_fluxes,1);
    #for flux_object in results_array

        # grab the flux -
        #flux_array = flux_object.flux_array;

        # cache -
        #flux_ensemble = [flux_ensemble flux_array];
    #end

    # cut off the zeros -
    #flux_ensemble = flux_ensemble[:,2:end];

    # compute the mean and std -
    #µ = mean(flux_ensemble, dims=2);
    #σ = (1/sqrt(number_of_samples))*std(flux_ensemble, dims=2);
    #flux_distribution = [µ σ];

    # return -
    return results_array
end

function maximize_specific_growth_rate_non_regulated(path_to_measurements_file::String;) #using calculate_optimal_flux_distribution

    # initalize -
    results_array = Array{VLOptimalFluxResult,1}()



    # ok, so lets sample ...


        # load the default data_dictionary -
        default_data_dictionary = generate_default_data_dictionary();

        # pass the default dictionary to a customization method -
        updated_data_dictionary = optimize_specific_growth_rate(default_data_dictionary);

        # update dictionary with experimental data?
        updated_data_dictionary = constrain_measured_fluxes(updated_data_dictionary, path_to_measurements_file);
        #updated_data_dictionary["flux_bounds_array"][83,:] .= 0.0

        # estimate the optimal flux distrubution -
        results_array = calculate_optimal_flux_distribution(updated_data_dictionary);

        flux_bounds_array = updated_data_dictionary["flux_bounds_array"]

        # build a return type -
        #fluxResult = VLOptimalFluxResult();
        #fluxResult.objective_value = objective_value;
        #fluxResult.flux_array = calculated_flux_array;
        #fluxResult.dual_array = dual_value_array;
        #fluxResult.uptake_array = uptake_array;
        #fluxResult.exit_flag = exit_flag;
        #fluxResult.status_flag = status_flag;
        #fluxResults.DRgenes = DRgenes;
        #fluxResults.constrainedRxns = constrainedRxns;
        #fluxResults.states = states;

        # get some problem setup information -
        #fluxResult.flux_bounds_array = regulated_data_dictionary["flux_bounds_array"];

        # cache -
        #push!(results_array, fluxResult);

    # how many flux are there?
    #number_of_fluxes = length(results_array[1].flux_array);

    # compute the flux values -
    #flux_ensemble = zeros(number_of_fluxes,1);
    #for flux_object in results_array

        # grab the flux -
        #flux_array = flux_object.flux_array;

        # cache -
        #flux_ensemble = [flux_ensemble flux_array];
    #end

    # cut off the zeros -
    #flux_ensemble = flux_ensemble[:,2:end];

    # compute the mean and std -
    #µ = mean(flux_ensemble, dims=2);
    #σ = (1/sqrt(number_of_samples))*std(flux_ensemble, dims=2);
    #flux_distribution = [µ σ];

    # return -
    return (results_array, flux_bounds_array)
end
