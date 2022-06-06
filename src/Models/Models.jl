module Models

using PolynomialRoots
using LinearAlgebra
using Gridap
using StaticArrays
using SparseArrays

# Export
export ShallowWater, ShallowWaterSolution
export solve, u₁, u₂
export ∂ₓu₁, ∂ₓu₂
export ∂ₓ²u₁, ∂ₓ²u₂
export ∂ₓ³u₁, ∂ₓ³u₂

export FiniteDepth, FiniteDepthSolution

export FiniteElementModel, FiniteDepthFEM

export ηₖ, ξₖ

export ϕₕ

# Import
import iceFEM.NonDimProblem: Ice, Fluid
import iceFEM.NonDimProblem: NonDimensionalProblem, non_dimensionalize, non_dimensionalize!, preallocate_matrices
import iceFEM.NonDimProblem: FreeBedrock, FreeClamped, FreeHinged, FreeFree

import iceFEM.DispersionEquations: dispersion_free_surface, dispersion_elastic_surface, solve_eigen_eb, dispersion_ice

include("ShallowWaterModel.jl")
include("FiniteDepthModel.jl")
include("FiniteElementModel.jl")

end
