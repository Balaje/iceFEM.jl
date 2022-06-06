using iceFEM
using Test

ice = Ice(922.5, 2e9, 0.33, 40000, 200);
fluid = Fluid(1025, 0, 9.8, 500, 40000);
ω = 2π/500
Aₚ = 9.8/(1im*ω)

partition = (400,40)
fe_model = FiniteElementModel(2, partition, 40, 8)
solFE = iceFEM.solve(ice, fluid, ω, FreeClamped(), fe_model; verbosity=0)

R₁ = solFE.linear_system[3][1]

@testset verbose= true "Checking the validity of finite element solution" begin
  @test abs(abs(R₁) - 1.) ≤ 1e-6
end

### Add tests for other beams here ...
