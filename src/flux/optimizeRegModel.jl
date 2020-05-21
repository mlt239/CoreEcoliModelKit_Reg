function optimizeRegModel(data_dictionary) #For now, we default to all of the states being true.
    rFBAsol1 = false.*ones(length(data_dictionary["regulatoryGenes"])); # initial state
    input1state = false.*ones(length(data_dictionary["regulatoryInputs1"]));
    input2state = false.*ones(length(data_dictionary["regulatoryInputs2"]));
    state1 = false.*ones(length(data_dictionary["regulatoryGenes"])+length(data_dictionary["regulatoryInputs1"])+length(data_dictionary["regulatoryInputs2"])); #vector for entire state, 0 and 1s
    optimizeRegModel2(data_dictionary, rFBAsol1, input1state, input2state, state1)
end

function optimizeRegModel(data_dictionary, initialRegState)
    if length(initialRegState) != (length(data_dictionary["regulatoryGenes"]) + length(data_dictionary["regulatoryInputs1"]) + length(data_dictionary["regulatoryInputs2"]))
        error("InitialRegState is invalid length")
    end
    state1 = initialRegState;
    rFBAsol1 = state1[1:length(data_dictionary["regulatoryGenes"])]
    input1state = state1[(length(data_dictionary["regulatoryGenes"])+1):(length(data_dictionary["regulatoryGenes"]) + length(data_dictionary["regulatoryInputs1"]))]
    input2state = state1[(length(data_dictionary["regulatoryGenes"]) + length(data_dictionary["regulatoryInputs1"]) + 1):(length(data_dictionary["regulatoryGenes"]) + length(data_dictionary["regulatoryInputs1"]) + length(data_dictionary["regulatoryInputs2"]))]
    optimizeRegModel2(data_dictionary, rFBAsol1, input1state, input2state, state1)
end

function optimizeRegModel2(data_dictionary, rFBAsol1, input1state, input2state, state1)
rFBAsols = [rFBAsol1];
input1states = [input1state]; #pay attention to the plural, this will hold all our states
input2states = [input2state];
states = [state1];

(rFBAsol2, final_input1state, final_input2state) = solveBooleanReg(data_dictionary, rFBAsol1, input1state, input2state)

state2 = [rFBAsol2;final_input1state;final_input2state]
push!(rFBAsols, rFBAsol2)
push!(input1states, final_input1state)
push!(input2states, final_input2state)
push!(states, state2)
println("Current gene iteration completed: ", length(states))
#now the looping begins
cycleReached = false
cycleStart = []
while ~cycleReached
    rFBAsol1 = rFBAsols[length(rFBAsols)]
    input1state = input1states[length(input1states)]
    input2state = input2states[length(input2states)]
    (rFBAsol2, final_input1state, final_input2state) = solveBooleanReg(data_dictionary, rFBAsol1, input1state, input2state)
    state2 = [rFBAsol2;final_input1state;final_input2state]
    push!(rFBAsols, rFBAsol2)
    push!(input1states, final_input1state)
    push!(input2states, final_input2state)
    push!(states, state2)
    println("Current gene iteration completed: ", length(states))
    for i = 1:(length(states)-1)
        if all(states[length(states)] == states[i])
            cycleStart = i
            println("Gene steady state reached after cycle: ", cycleStart)
        end
    end
    if ~isempty(cycleStart)
        cycleReached = true
    end
end

#now lets do the final FBAs for the steady state regulatory network
FBAsols = []
DRgenes = []
constrainedRxns = []
objsols = []
regulated_data_dicts = []
dv_arrays = []
up_arrays = []
exit_flags = []
status_flags = []
k = 0
for i = cycleStart:(length(states)-1)
    k = k+1
    dgenes = []
    for j = 1:length(data_dictionary["regulatoryGenes"])
        if rFBAsols[i][j] == false
            push!(dgenes, data_dictionary["regulatoryGenes"][j])
        end
    end
    dgenes = intersect(data_dictionary["list_of_gene_symbols"], dgenes)
    (regulated_data_dictionary, hasEffect, crxns) = deleteModelGenes(data_dictionary, dgenes)
    #regulated_data_dictionary = optimize_specific_growth_rate(regulated_data_dictionary) #placeholder
    (objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag) = calculate_optimal_flux_distribution(regulated_data_dictionary); #change this to 2 version if standalone
    #println("Look at this:", crxns)
    push!(FBAsols, calculated_flux_array)
    push!(DRgenes, dgenes)
    push!(constrainedRxns, crxns)
    push!(objsols, objective_value)
    push!(regulated_data_dicts, regulated_data_dictionary)
    push!(dv_arrays, dual_value_array)
    push!(up_arrays, uptake_array)
    push!(exit_flags, exit_flag)
    push!(status_flags, status_flag)


#return(regulated_data_dictionary, objsols, FBAsols, dual_value_array, uptake_array, DRgenes, constrainedRxns, states, exit_flag, status_flag)
end
return(regulated_data_dicts, objsols, FBAsols, dv_arrays, up_arrays, DRgenes, constrainedRxns, states, exit_flags, status_flags)
end
