using iceFEM
using Plots

L = 10000
ρᵢ = 922.5
Eᵢ = 2e9
ν = 0.33
h = 200
ice = Ice(ρᵢ, Eᵢ, ν, L, h)
ρₒ = 1025.0
g = 9.8
H = 500
k₀ = 1e6
ϵ = 0.7
fluid = Fluid(ρₒ, k₀, g, H, ϵ*L)

ω = 2π/100
# Solve Bedrock problem using FEM
fe_model = FiniteElementModel(2, (100,20), 10, 5) # dim, partition, nev, NModes
sol_bedrock = iceFEM.solve(ice, fluid, ω, FreeBedrock(), fe_model)

# Change shelf length and run for clamped BC using Finite Depth
ice = Ice(ρᵢ, Eᵢ, ν, ϵ*L, h)
fluid = Fluid(ρₒ, 0, g, H, ϵ*L)
sol_clamped_emm = iceFEM.solve(ice, fluid, ω, FreeClamped(), FiniteDepth(5))

# Solve the same problem using FEM with clamped BC
# NOTE: Switch on verbose mode.
sol_clamped = iceFEM.solve(ice, fluid, ω, FreeClamped(), fe_model, verbosity = 1)

# Plot all the solutions
x₁ = LinRange(0, sol_bedrock.ndp.geo[4], 200)
x₂ = LinRange(sol_bedrock.ndp.geo[4], sol_bedrock.ndp.geo[1], 200)
plt = plot(x₁, abs.(u₁(x₁, sol_bedrock)), label="\$x < x_g\$ FEM solution", line=(:solid, 2))
plot!(plt, x₂, abs.(u₂(x₂, sol_bedrock)), label="\$x > x_g\$ FEM solution", line=(:solid, 2))
plot!(plt, x₁, abs.(u₁(x₁, sol_clamped_emm)), label="Clamped at \$x=x_g\$ EMM solution", line=(:dash, 1.5))
plot!(plt, x₁, abs.(u₁(x₁, sol_clamped)), label="Clamped at \$x=x_g\$ FEM solution", line=(:dash, 1.5))

xlabel!(plt, "\$\\bar{x}\$ (in non-dim)")
ylabel!(plt, "\$|u|\$ (in m)")
