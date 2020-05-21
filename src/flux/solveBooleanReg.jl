function solveBooleanReg(data_dictionary, initial_state, initial_input1, initial_input2)

#bring in the data that I need to evaluate the regulatory states
regulatoryRules = data_dictionary["regulatoryRules"]
regulatoryInputs1 = data_dictionary["regulatoryInputs1"]
regulatoryInputs2 = data_dictionary["regulatoryInputs2"]
flux_bounds = data_dictionary["flux_bounds_array"]
regulatoryGenes = data_dictionary["regulatoryGenes"]
genes = data_dictionary["list_of_gene_symbols"]
flux_names = data_dictionary["list_of_reaction_name_strings"]
final_input1 = false.*ones(length(data_dictionary["regulatoryInputs1"])); #initialize
s = data_dictionary["stoichiometric_matrix"]

initial_state = initial_state .== 1
initial_input1 = initial_input1 .== 1
initial_input2 = initial_input2 .== 1
initial_regstate = vcat(initial_state, initial_input1, initial_input2) #Full list of states

#This determines external metabolite levels using flux bounds (input 1)
for i = 1:length(regulatoryInputs1)
    flux_n = regulatoryInputs1[i]; #get the metabolite in the first regulatory input
    metID = findfirst(isequal(flux_n), data_dictionary["list_of_metabolite_symbols"]) #find the index of that metabolite
    potfluxID = findall(isequal(~0), s[metID[1],:]) #find the fluxes associated with that metabolite
#so now we only need to select actual exchanges, meaning only one -1 or 1.
    for j = 1:length(potfluxID)
        if abs(sum(s[:,potfluxID[j]])) == 1
            rxnID = potfluxID[j]
            #println(rxnID)
                if flux_bounds[rxnID, 1] < 0 #if exchange is allowed
                    final_input1[i] = true;
                    #println(flux_bounds[rxnID,1])
                    #println("True")
                else
                    final_input1[i] = false; #if exchange is not allowed
                    #println(flux_bounds[rxnID,1])
                    #println("false")
                end
        end
    end

end #now we have new input1 states

#apply initial regulatory state to model and obtain FBA
drGenes = []
for i = 1:length(regulatoryGenes)
    if initial_state[i] == false
        push!(drGenes, regulatoryGenes[i])
    end
end
drGenes = intersect(genes,drGenes)
data_dictionary_deleted = deleteModelGenes(data_dictionary, drGenes)[1]
#data_dictionary_deleted = optimize_specific_growth_rate(data_dictionary_deleted) #placeholder for now
(objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag) = calculate_optimal_flux_distribution(data_dictionary_deleted); #change this to 2 version if standalone
#println("Succ:", calculated_flux_array[87])
#Now, we need to check regulatory inputs that are dependent on fluxes
final_input2 = true.*ones(length(data_dictionary["regulatoryInputs2"]))
#println(calculated_flux_array)
fluxInd = indexin(regulatoryInputs2, flux_names) #figure out where these inputs are in the list of fluxes
for i = 1:length(regulatoryInputs2)
    if calculated_flux_array[fluxInd[i]] == 0 #if there was no flux in the regulatory input 2 flux
        final_input2[i] = false
    end
end

#Now that we have regulatory inputs 1 and 2 set, it's time to set the state of genes for next time

geneStates = calculate_full_rules_vector(data_dictionary, initial_regstate)
#geneStates = geneStates .== 1
final_state = geneStates #the final state of the regulatory system

return(final_state,final_input1,final_input2)
#so the next time around, the final_state will differ from the initial state as
#we now are inputting new inputs1 and inputs2.
end
