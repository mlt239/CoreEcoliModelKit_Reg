#This function will actually constrain the model based on deleted genes

function deleteModelGenes(data_dictionary, drGenes)
genes = deepcopy(data_dictionary["list_of_gene_symbols"]) #need to update this later
genes_minus_deletions = [] #our genelist will eventually be marked
rxnGeneMat = data_dictionary["rxnGeneMat"] #need this for later
flux_bounds_array = deepcopy(data_dictionary["flux_bounds_array"]) #will need to constrain certain fluxes
rules = deepcopy(data_dictionary["irules"])
hasEffect = false #flag to let us know if something actually changes
constrainRxnNames = []


for i = 1:length(genes)
    push!(genes_minus_deletions, replace(genes[i], "_deleted" => ""))
end

#Check if there's genes that don't exist in the deleted genes list
geneInd = indexin(drGenes, genes_minus_deletions)
flag = false
for i = 1:length(geneInd)
    if geneInd[i] == nothing
        flag = true
    end
end
#println(flag)

if ~flag #if our list of deleted genes is in the model, proceed
deletedGeneNames = []
    for i in geneInd
        append!(deletedGeneNames, genes[i])
        genes[i] = genes_minus_deletions[i]*"_deleted" #marks the genes for deletion
    end

    #lets grab our list of rxns that are associated with our deleted genes
    rxnInd = []
    for i in geneInd
        append!(rxnInd, findall(rxnGeneMat[:,i] .== 1))
    end # now rxnInd inludes rxns that are affected by the deleted genes

    if ~isempty(rxnInd) #if affected rxns
        x = true.*ones(length(genes)) #represent genes as booleans
        for i in geneInd #set all of the deleted genes to false for the boolean
            x[i] = false
        end

        x = x .== 1 #convert x to a bit array that can be evaluated as a boolean
        constrainRxns = false.*ones(length(rxnInd))
        evaluatedRules = calculate_gene_rules_vector(data_dictionary, x)
        evaluatedRules = evaluatedRules .== 1
        #now it's time to figure out which rxns are constrained
        for j = 1:length(rxnInd)
            if ~isempty(evaluatedRules[rxnInd[j]]) #if the rules are empty, theres no regulation
                if ~evaluatedRules[rxnInd[j]] #if the rule eval turns out false..
                    constrainRxns[j] = true #this rxn will be constrained
                end
            end
        end
        if any(constrainRxns .== true) #if we have constrainedRxns
            for k = 1:length(constrainRxns)
                if constrainRxns[k] == true #if the particular index is true
                    push!(constrainRxnNames, data_dictionary["list_of_reaction_name_strings"][rxnInd[k]])
                    flux_bounds_array[rxnInd[k],:] .= 0 #silence by setting flux to 0
                end
            end
            hasEffect = true
        end
    end

    #output updated gene list, flux bounds,
    deleted_data_dictionary = deepcopy(data_dictionary)
    deleted_data_dictionary["list_of_gene_symbols"] = genes
    deleted_data_dictionary["flux_bounds_array"] = flux_bounds_array
    trueconstrainRxnNames = unique(constrainRxnNames)
    return (deleted_data_dictionary, hasEffect, trueconstrainRxnNames, deletedGeneNames)

else
    error("a gene is not in the model!")
end
end


#function buildRxnGeneMat(data_dictionary)
#rxnGeneMat = false.*ones(data_dictionary["number_of_reactions"], length(data_dictionary["list_of_gene_symbols"])

#end
