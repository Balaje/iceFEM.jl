module DispersionEquations

using Polynomials
using Roots
using LinearAlgebra

# Import
import iceFEM.NonDimProblem: Ice, Fluid
import iceFEM.NonDimProblem: NonDimensionalProblem, non_dimensionalize
import iceFEM.NonDimProblem: FreeBedrock, FreeClamped, FreeHinged

# Export
export dispersion_elastic_surface, dispersion_free_surface, solve_eigen_eb

include("DispersionFreeSurface.jl")
include("DispersionElasticSurface.jl")
include("BeamDispersionEquations.jl")

end
