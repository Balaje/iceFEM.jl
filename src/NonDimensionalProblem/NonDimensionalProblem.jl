module NonDimProblem

# Export
export Ice, Fluid, NonDimensionalProblem
export FreeBedrock, FreeClamped, FreeHinged, FreeFree

export non_dimensionalize, non_dimensionalize!, preallocate_matrices

include("Input.jl")
end
