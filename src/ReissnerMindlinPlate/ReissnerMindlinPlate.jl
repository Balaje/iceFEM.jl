module ReissnerMindlinPlate

using PolynomialRoots
using LinearAlgebra
using StaticArrays
using NLsolve

# Export
export ReissnerMindlinIce
export non_dimensionalize
export dispersion_ice
export solve

# Import
import iceFEM.NonDimProblem: Fluid, Ice
import iceFEM.NonDimProblem: NonDimensionalProblem, non_dimensionalize
import iceFEM.NonDimProblem: FreeBedrock, FreeClamped, FreeHinged, FreeFree

import iceFEM.DispersionEquations: dispersion_free_surface, dispersion_ice

import iceFEM.Models: FiniteDepth, solve, FiniteDepthSolution


struct ReissnerMindlinIce <: Any end


include("NonDimensionalProblem.jl");
include("DispersionEquation.jl");
include("FiniteDepthModel.jl");

end
