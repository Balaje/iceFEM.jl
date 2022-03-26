module DispersionEquations

# Export
export dispersion_elastic_surface, dispersion_free_surface

include("DispersionFreeSurface.jl")
include("DispersionElasticSurface.jl")

end
