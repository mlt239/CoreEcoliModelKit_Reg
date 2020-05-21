# declare a return type for the flux variability calculation -
mutable struct VLFluxVariabilityResult

    optimal_flux_array::Array{Float64,1}
    flux_variability_array::Array{Float64,1}

    # constructor -
    function VLFluxVariabilityResult()
        this = new()
    end
end

# declare a return type -
mutable struct VLOptimalFluxResult

    objective_value::Float64
    flux_array::Array{Float64,1}
    dual_array::Array{Float64,1}
    uptake_array::Array{Float64,1}
    exit_flag::Int64
    status_flag::Int64
    flux_bounds_array::Array{Float64,2}

    # constructor -
    function VLOptimalFluxResult()
        this = new()
    end
end

mutable struct Gene

    # information about a Gene -
    symbol::AbstractString
    organism_code::AbstractString


    function Gene()
        this = new()
    end
end

mutable struct Promoter

    # information about a Gene -
    symbol::AbstractString
    organism_code::AbstractString


    function Promoter()
        this = new()
    end
end
