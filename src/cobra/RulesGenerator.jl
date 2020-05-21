function preParse(data_dictionary)
    regulatoryRules = deepcopy(data_dictionary["regulatoryRules"])
    regulatoryFull = vcat(data_dictionary["regulatoryGenes"], data_dictionary["regulatoryInputs1"], data_dictionary["regulatoryInputs2"])
    rulesList = Array{Any, 2}(nothing, size(regulatoryRules))
    for k = 1:length(regulatoryRules)
        fields = split(regulatoryRules[k], r"([\s~|&()])", keepempty = false)
        newfields = copy(fields)
        for j = 1:length(fields)
            word = fields[j]
            if word == "true"
                newfields[j] = ""
            elseif word == "false"
                newfield[j] = "false"
            else
                if ~isempty(findall(regulatoryFull .== word))
                    index = findall(regulatoryFull .== word)
                    newfields[j] = "initial_regstate("*string(index[1][1])*")"
                else
                    newfields[j] = ""
                end
            end
        end
            newRule = regulatoryRules[k]
            for j = 1:length(fields)
                newRule = replace(newRule, fields[j] => newfields[j])
            end
            rulesList[k] = newRule
    end
    for i = 1:length(rulesList)
        rulesList[i] = replace(rulesList[i], "||" => "|")
        rulesList[i] = replace(rulesList[i], "&&" => "&")
    end

    return rulesList
end

function process(expression::Expr)

    # new expression -
    new_expression = copy(expression)

    # get the args for this expression -
    args_array = expression.args
    number_of_arguments = length(args_array)
    for arg_index = 1:number_of_arguments

        # get the arg -
        arg = args_array[arg_index]

        if typeof(arg) == Symbol
            if arg == :|
                new_expression.args[arg_index] = :max
            elseif arg == :&
                new_expression.args[arg_index] = :min
            elseif arg == :x
                new_expression.args[arg_index] = :fn
            elseif arg == :initial_regstate
                new_expression.args[arg_index] = :fn
            elseif arg == :~
                new_expression.args[arg_index] = :~
            end
        elseif typeof(arg) == Expr
            new_expression.args[arg_index] = process(arg)
        end
    end

    return new_expression
end

function generate_rules_function(cobra_dictionary, path_to_output_rules_function::String)

    # Initialize -
    expression_array = Array{Any,1}()

    # Get the rules -
    rules_array = cobra_dictionary["rules"]

    # how many rules do we have?
    number_of_rules = length(rules_array)

    # create an array of expressions -
    for rule_string in rules_array

        # build an expression -
        local_expr = Meta.parse(rule_string)

        # process this expression -
        if local_expr == nothing
            push!(expression_array, nothing)
        else

            # work on this tree -
            converted_expression = process(local_expr)

            # grab -
            push!(expression_array, converted_expression)
        end
    end



    # dump reaction array to disk -
    open(path_to_output_rules_function, "w") do f

        # buffer -
        buffer = String[]
        tmp_value = "function calculate_rules_vector(data_dictionary::Dict{String,Any}, value_array::Array{Float64,1})"
        push!(buffer, tmp_value)
        push!(buffer, "")
        push!(buffer, "\t# Initialize - ")
        push!(buffer, "\tv = zeros($(number_of_rules))")
        push!(buffer, "")

        # write the buffer -
        for line_item in buffer
            write(f,"$(line_item)\n")
        end

        # define some *internal* functions -
        default_fn_buffer = String[]
        push!(default_fn_buffer,"\tfunction estimate_default(data_dictionary)")
        push!(default_fn_buffer,"\t\treturn 1.0")
        push!(default_fn_buffer,"\tend")
        push!(default_fn_buffer, "\tdefault() = estimate_default(data_dictionary)")
        push!(default_fn_buffer,"")

        # write the default_fn_buffer -
        for line_item in default_fn_buffer
            write(f,"$(line_item)\n")
        end

        # fn buffer -
        fn_buffer = String[]
        push!(fn_buffer, "\tfunction calculate_fn_value(x, value_array)")
        push!(fn_buffer, "\t\treturn value_array[x]")
        push!(fn_buffer,"\tend")
        push!(fn_buffer, "\tfn(x) = calculate_fn_value(x, value_array)")
        push!(fn_buffer,"")

        # write the default_fn_buffer -
        for line_item in fn_buffer
            write(f,"$(line_item)\n")
        end

        # write the rules -
        counter = 1
        for rule in expression_array

            line_item = ""

            if rule == nothing
                line_item = "\tv[$(counter)] = default()"
            else

                # get the rule string -
                rule_string = repr(rule)
                line_item = "\tv[$(counter)] = $(rule_string[3:end-1])"
            end

            counter = counter + 1
            write(f,"$(line_item)\n")
        end

        # return -
        write(f,"\n")
        write(f,"\treturn v\n")

        # write the end -
        write(f, "end\n")
    end
end

function generate_regulatory_function(cobra_dictionary, path_to_output_rules_function::String)

    # Initialize -
    expression_array = Array{Any,1}()

    # Get the rules -
    full_rules_array = preParse(cobra_dictionary)

    # how many rules do we have?
    number_of_rules = length(full_rules_array)

    # create an array of expressions -
    for rule_string in full_rules_array

        # build an expression -
        local_expr = Meta.parse(rule_string)
        println(local_expr)
        # process this expression -
        if local_expr == nothing
            push!(expression_array, nothing)
        else

            # work on this tree -
            converted_expression = process(local_expr)

            # grab -
            push!(expression_array, converted_expression)
        end
    end



    # dump reaction array to disk -
    open(path_to_output_rules_function, "w") do f

        # buffer -
        buffer = String[]
        tmp_value = "function calculate_full_rules_vector(data_dictionary::Dict{String,Any}, value_array::BitArray{2})"
        push!(buffer, tmp_value)
        push!(buffer, "")
        push!(buffer, "\t# Initialize - ")
        push!(buffer, "\tv = zeros($(number_of_rules))")
        push!(buffer, "")

        # write the buffer -
        for line_item in buffer
            write(f,"$(line_item)\n")
        end

        # define some *internal* functions -
        default_fn_buffer = String[]
        push!(default_fn_buffer,"\tfunction estimate_default(data_dictionary)")
        push!(default_fn_buffer,"\t\treturn 1.0")
        push!(default_fn_buffer,"\tend")
        push!(default_fn_buffer, "\tdefault() = true")
        push!(default_fn_buffer,"")

        # write the default_fn_buffer -
        for line_item in default_fn_buffer
            write(f,"$(line_item)\n")
        end

        # fn buffer -
        fn_buffer = String[]
        push!(fn_buffer, "\tfunction calculate_fn_value(x, value_array)")
        push!(fn_buffer, "\t\treturn value_array[x]")
        push!(fn_buffer,"\tend")
        push!(fn_buffer, "\tfn(x) = calculate_fn_value(x, value_array)")
        push!(fn_buffer,"")

        # write the default_fn_buffer -
        for line_item in fn_buffer
            write(f,"$(line_item)\n")
        end

        # write the rules -
        counter = 1
        for rule in expression_array

            line_item = ""

            if rule == nothing
                line_item = "\tv[$(counter)] = default()"
            else

                # get the rule string -
                rule_string = repr(rule)
                line_item = "\tv[$(counter)] = $(rule_string[3:end-1])"
            end

            counter = counter + 1
            write(f,"$(line_item)\n")
        end

        # return -
        write(f,"\n")
        write(f,"\treturn v\n")

        # write the end -
        write(f, "end\n")
    end
end
