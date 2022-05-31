module ComplexResonances

using PolynomialRoots
using LinearAlgebra
using Interpolations
using GenericLinearAlgebra
using LinearAlgebra

# export modules
export InterpolateFreqDomain, FreqSpace
export computeResonanceFrequency

# import modules
import iceFEM.NonDimProblem: Fluid, Ice
import iceFEM.NonDimProblem: NonDimensionalProblem, non_dimensionalize
import iceFEM.NonDimProblem: FreeBedrock, FreeClamped, FreeHinged, FreeFree

import iceFEM.Models: FiniteDepth, solve, FiniteDepthSolution

include("InterpolateFreq.jl")
include("ResonanceFreq.jl")

end
