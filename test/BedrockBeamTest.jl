using iceFEM
using Test

#######################
# Create new problem
#######################
ice = Ice(922.5, 2e9, 0.33, 40000, 200);
fluid = Fluid(1025, 1e6, 9.8, 500, 0.7*40000);
ω = 2π/500
Aₚ = 9.8/(1im*ω)

#######################
# Solve bedrock beam under
# 1) Shallow-water conditions
# 2) Finite depth conditions
#######################
sol_sw = solve(ice, fluid, ω, FreeBedrock(), ShallowWater())
sol_fd = solve(ice, fluid, ω, FreeBedrock(), FiniteDepth())
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

xg = sol_sw.ndp.geo[4]
@testset verbose = true "Checking whether the bedrock coupling conditions are satisfied at the Grounding line x=xg" begin
  @testset "Shallow water..." begin
    @test u₁(xg, sol_sw) ≈ u₂(xg, sol_sw)
    @test ∂ₓu₁(xg, sol_sw) ≈ ∂ₓu₂(xg, sol_sw)
    @test ∂ₓ²u₁(xg, sol_sw) ≈ ∂ₓ²u₂(xg, sol_sw)
    @test ∂ₓ³u₁(xg, sol_sw) ≈ ∂ₓ³u₂(xg, sol_sw)
  end
  @testset "Finite Depth..." begin
    @test u₁(xg, sol_fd) ≈ u₂(xg, sol_fd)
    @test ∂ₓu₁(xg, sol_fd) ≈ ∂ₓu₂(xg, sol_fd)
    @test ∂ₓ²u₁(xg, sol_fd) ≈ ∂ₓ²u₂(xg, sol_fd)
    @test ∂ₓ³u₁(xg, sol_fd) ≈ ∂ₓ³u₂(xg, sol_fd)
  end
end

@testset verbose=true "Checking the non-local boundary condition at the Grounding line" begin
  @testset "Shallow Water..." begin
    𝐴 = (-sol_sw.p[1])*(-sol_sw.p[2])
    𝐵 = (-sol_sw.p[1] + -sol_sw.p[2])
    𝐶 = ( (-sol_sw.p[1])^2 + (-sol_sw.p[1])*(-sol_sw.p[2]) + (-sol_sw.p[2])^2 )
    @test -𝐴*(u₁(xg, sol_sw)) + (𝐵 )*(∂ₓu₁(xg,sol_sw)) ≈ ∂ₓ²u₁(xg, sol_sw)
    @test (-𝐴*𝐵)*(u₁(xg, sol_sw)) + (𝐶 )*(∂ₓu₁(xg,sol_sw)) ≈ ∂ₓ³u₁(xg, sol_sw)
  end
  @testset "Finite Depth..." begin
  𝐴 = (-sol_fd.p[1])*(-sol_fd.p[2])
  𝐵 = (-sol_fd.p[1] + -sol_fd.p[2])
  𝐶 = ( (-sol_fd.p[1])^2 + (-sol_fd.p[1])*(-sol_fd.p[2]) + (-sol_fd.p[2])^2 )
    @test -𝐴*(u₁(xg, sol_fd)) + (𝐵 )*(∂ₓu₁(xg,sol_fd)) ≈ ∂ₓ²u₁(xg, sol_fd)
    @test (-𝐴*𝐵)*(u₁(xg, sol_fd)) + (𝐶 )*(∂ₓu₁(xg,sol_fd)) ≈ ∂ₓ³u₁(xg, sol_fd)
  end
end
