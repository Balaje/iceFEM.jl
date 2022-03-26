module Models

using Polynomials

# Export
export ShallowWater, ShallowWaterSolution
export solve, u₁, u₂
export ∂ₓu₁, ∂ₓu₂
export ∂ₓ²u₁, ∂ₓ²u₂
export ∂ₓ³u₁, ∂ₓ³u₂

export FiniteDepth, FiniteDepthSolution

export ηₖ, ξₖ

# Import
import iceFEM.NonDimProblem: Ice, Fluid
import iceFEM.NonDimProblem: NonDimensionalProblem, non_dimensionalize
import iceFEM.NonDimProblem: FreeBedrock, FreeClamped, FreeHinged

import iceFEM.DispersionEquations: dispersion_free_surface, dispersion_elastic_surface, solve_eigen_eb

include("ShallowWaterModel.jl")
include("FiniteDepthModel.jl")
include("FiniteElementModel.jl")

end
