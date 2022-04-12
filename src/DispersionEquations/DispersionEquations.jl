module DispersionEquations

using PolynomialRoots
using LinearAlgebra

# Import
import iceFEM.NonDimProblem: Ice, Fluid
import iceFEM.NonDimProblem: NonDimensionalProblem, non_dimensionalize
import iceFEM.NonDimProblem: FreeBedrock, FreeClamped, FreeHinged

# Export
export dispersion_elastic_surface, dispersion_free_surface, solve_eigen_eb
export dispersion_ice

include("DispersionFreeSurface.jl")
include("DispersionElasticSurface.jl")
include("BeamDispersionEquations.jl")
include("DispersionIce.jl")

end
