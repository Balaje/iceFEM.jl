module iceFEM

using Polynomials
using Plots
using LaTeXStrings
using LinearAlgebra
using Roots

include("./NonDimensionalProblem/NonDimensionalProblem.jl")
include("./DispersionEquations/DispersionEquations.jl")
include("./Models/Models.jl")

# Export
using iceFEM.NonDimProblem: Ice, Fluid, NonDimensionalProblem
using iceFEM.NonDimProblem: FreeBedrock, FreeClamped, FreeHinged
using iceFEM.NonDimProblem: non_dimensionalize

using iceFEM.DispersionEquations: dispersion_elastic_surface, dispersion_free_surface
using iceFEM.DispersionEquations: solve_eigen_eb

using iceFEM.Models: ShallowWater, ShallowWaterSolution, solve
using iceFEM.Models: u₁, ∂ₓu₁, ∂ₓ²u₁, ∂ₓ³u₁
using iceFEM.Models: u₂, ∂ₓu₂, ∂ₓ²u₂, ∂ₓ³u₂
using iceFEM.Models: FiniteDepth, FiniteDepthSolution
using iceFEM.Models: FiniteDepthFEM, FiniteElementModel

export Ice, Fluid, NonDimensionalProblem
export FreeBedrock, FreeClamped, FreeHinged
export non_dimensionalize
export dispersion_free_surface, dispersion_elastic_surface, solve_eigen_eb
export ShallowWater, ShallowWaterSolution, solve
export FiniteDepth, FiniteDepthSolution
export FiniteDepthFEM, FiniteElementModel
export u₁, ∂ₓu₁, ∂ₓ²u₁, ∂ₓ³u₁
export u₂, ∂ₓu₂, ∂ₓ²u₂, ∂ₓ³u₂

end # module
