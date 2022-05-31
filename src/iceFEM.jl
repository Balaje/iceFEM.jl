module iceFEM

using PolynomialRoots
using Plots
using LaTeXStrings
using LinearAlgebra
using Roots

include("./NonDimensionalProblem/NonDimensionalProblem.jl")
include("./DispersionEquations/DispersionEquations.jl")
include("./Models/Models.jl")
include("./ReissnerMindlinPlate/ReissnerMindlinPlate.jl")
include("./ComplexResonances/ComplexResonances.jl")

# Export
using iceFEM.NonDimProblem: Ice, Fluid, NonDimensionalProblem
using iceFEM.NonDimProblem: FreeBedrock, FreeClamped, FreeHinged, FreeFree
using iceFEM.NonDimProblem: non_dimensionalize

using iceFEM.DispersionEquations: dispersion_elastic_surface, dispersion_free_surface
using iceFEM.DispersionEquations: dispersion_ice
using iceFEM.DispersionEquations: solve_eigen_eb

using iceFEM.Models: ShallowWater, ShallowWaterSolution, solve
using iceFEM.Models: u₁, ∂ₓu₁, ∂ₓ²u₁, ∂ₓ³u₁
using iceFEM.Models: u₂, ∂ₓu₂, ∂ₓ²u₂, ∂ₓ³u₂
using iceFEM.Models: FiniteDepth, FiniteDepthSolution
using iceFEM.Models: FiniteDepthFEM, FiniteElementModel

using iceFEM.ReissnerMindlinPlate: ReissnerMindlinIce
using iceFEM.ReissnerMindlinPlate: non_dimensionalize, dispersion_ice

using iceFEM.ComplexResonances: InterpolateFreqDomain, FreqSpace
using iceFEM.ComplexResonances: computeResonanceFrequency

export Ice, Fluid, NonDimensionalProblem
export FreeBedrock, FreeClamped, FreeHinged, FreeFree
export non_dimensionalize
export dispersion_free_surface, dispersion_elastic_surface, solve_eigen_eb
export dispersion_ice
export ShallowWater, ShallowWaterSolution, solve
export FiniteDepth, FiniteDepthSolution
export FiniteDepthFEM, FiniteElementModel
export u₁, ∂ₓu₁, ∂ₓ²u₁, ∂ₓ³u₁
export u₂, ∂ₓu₂, ∂ₓ²u₂, ∂ₓ³u₂

export ReissnerMindlinIce

export InterpolateFreqDomain, FreqSpace, computeResonanceFrequency

end # module
