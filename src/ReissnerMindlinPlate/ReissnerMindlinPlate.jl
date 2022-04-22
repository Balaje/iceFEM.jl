module ReissnerMindlinPlate

using PolynomialRoots
using LinearAlgebra
using StaticArrays

# Export
export ReissnerMindlinIce
export non_dimensionalize

# Import
import iceFEM.NonDimProblem: Fluid, Ice
import iceFEM.NonDimProblem: NonDimensionalProblem, non_dimensionalize
import iceFEM.NonDimProblem: FreeBedrock, FreeClamped, FreeHinged, FreeFree

import iceFEM.DispersionEquations: dispersion_free_surface

import iceFEM.Models: FiniteDepth


struct ReissnerMindlinIce <: Any end


include("NonDimensionalproblem.jl");

end
