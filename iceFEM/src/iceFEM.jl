module iceFEM

import Polynomials
import Plots
import LaTeXStrings
import LinearAlgebra
import Roots


include("./NonDimensionalProblem/NonDimensionalProblem.jl")
include("./DispersionEquations/DispersionEquations.jl")
include("./Models/Models.jl")

# Export
using iceFEM.NonDimProblem: Ice, Fluid, NonDimensionalProblem
using iceFEM.NonDimProblem: FreeBedrock, FreeClamped, FreeHinged
using iceFEM.NonDimProblem: non_dimensionalize

using iceFEM.DispersionEquations: dispersion_elastic_surface, dispersion_free_surface

using iceFEM.Models: ShallowWater, ShallowWaterSolution, solve
using iceFEM.Models: FiniteDepth, FiniteDepthSolution


export Ice, Fluid, NonDimensionalProblem
export FreeBedrock, FreeClamped, FreeHinged
export non_dimensionalize
export dispersion_free_surface, dispersion_elastic_surface
export ShallowWater, ShallowWaterSolution, solve
export FiniteDepth, FiniteDepthSolution

end # module
