using iceFEM
using Test

#######################
# Create new problem
#######################
ice = Ice(922.5, 2e9, 0.33, 40000, 200);
fluid = Fluid(1025, 0, 9.8, 500, 40000);
ω = 2π/500
Aₚ = 9.8/(1im*ω)

#######################
# Solve clamped beam under
# 1) Shallow-water conditions
# 2) Finite depth conditions
#######################
sol_sw = solve(ice, fluid, ω, FreeClamped(), ShallowWater())
sol_fd = solve(ice, fluid, ω, FreeClamped(), FiniteDepth())
@testset verbose = true "Checking whether energy conservation is verified" begin
  @testset "Shallow water..." begin @test abs(sol_sw.a₀[1]/Aₚ) ≈ 1 end
  @testset "Finite Depth..." begin @test abs(sol_fd.aₘ[1]/Aₚ) ≈ 1 end
end


@testset verbose = true "Checking whether the free conditions are satisfied at x=0" begin
  @testset "Shallow water..." begin
    @test abs.(∂ₓ²u₁(0, sol_sw))[1] ≤ 1e-12
    @test abs.(∂ₓ³u₁(0, sol_sw))[1] ≤ 1e-12
  end
  @testset "Finite Depth..." begin
    @test abs.(∂ₓ²u₁(0, sol_fd))[1] ≤ 1e-12
    @test abs.(∂ₓ³u₁(0, sol_fd))[1] ≤ 1e-12
  end
end

LL = sol_sw.ndp.geo[1]
@testset verbose = true "Checking whether the clamped conditions are satisfied at x=L" begin
  @testset "Shallow water..." begin
    @test abs.(u₁(LL, sol_sw))[1] ≤ 1e-12
    @test abs.(∂ₓu₁(LL, sol_sw))[1] ≤ 1e-12
  end
  @testset "Finite Depth..." begin
    @test abs.(u₁(LL, sol_fd))[1] ≤ 1e-12
    @test abs.(∂ₓu₁(LL, sol_fd))[1] ≤ 1e-12
  end
end
