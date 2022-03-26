include("ShallowWaterModel.jl")

##################################################
# Ice-shelf and fluid/bedrock parameters
##################################################
L = 40000
ρᵢ = 922.5
Eᵢ = 2e9
ν = 0.33
h = 200
ice = Ice(ρᵢ, Eᵢ, ν, L, h)

ρₒ = 1025.0
g = 9.8
H = 500
fluid = Fluid(ρₒ, 0, g, H, 0) # Set the bedrock parameters = 0

ω = 2π/8000
#################################################
# Eg 1: Solve one frequency domain problem for clamped ice-shelf
#################################################
BeamType = FreeClamped()
WaterType = ShallowWater()
sol = solve(ice, fluid, ω, BeamType, WaterType)
Aₚ = g/(1im*ω)
@show abs(sol.a₀[1]/Aₚ)
# Visualize the solution
LL = sol.ndp.geo[1]
x₁ = 0:0.05:LL
U₁ = u₁(x₁, sol)
p1 = plot(x₁, abs.(U₁), legend=false)
xlabel!(p1, "\$x\$ (Non-Dim)")
ylabel!(p1, "\$|u|\$ (in m)")
title!(p1, "Clamped \$x = L\$, Wave Period = \$"*string(round(2π/ω, digits=4))*"\$ s")

#################################################
# Eg 1: Solve one frequency domain problem for hinged ice-shelf
#################################################
BeamType = FreeHinged()
sol = solve(ice, fluid, ω, BeamType, WaterType)
Aₚ = g/(1im*ω)
@show abs(sol.a₀[1]/Aₚ)
# Visualize the solution
LL = sol.ndp.geo[1]
x₁ = 0:0.05:LL
U₁ = u₁(x₁, sol)
p2 = plot(x₁, abs.(U₁), legend=false)
xlabel!(p2, "\$x\$ (Non-Dim)")
ylabel!(p2, "\$|u|\$ (in m)")
title!(p2, "Hinged \$x = L\$")

#################################################
# Eg 3: Solve one frequency domain problem for bedrock ice-shelf
# Note: Now the length should change along with the bedrock parameters
#################################################
L = 40000/0.7
ice = Ice(ρᵢ, Eᵢ, ν, L, h)
k₀ = 1e9 # Stiffness of the bedrock
x₀ = 0.7*L # For the grounding line problem
fluid = Fluid(ρₒ, k₀, g, H, x₀)
BeamType = FreeBedrock()
sol = solve(ice, fluid, ω, BeamType, WaterType)
Aₚ = g/(1im*ω)
@show abs(sol.a₀[1]/Aₚ)
# Visualize the solution
xg = sol.ndp.geo[4]
LL = sol.ndp.geo[1]
x₁ = 0:0.01:xg
x₂ = xg:0.01:LL
U₁ = u₁(x₁, sol)
U₂ = u₂(x₂, sol)
p3 = plot(x₁, abs.(U₁), label="\$x < x_g\$")
plot!(p3, x₂, abs.(U₂), label="\$ x > x_g\$")
xlabel!(p3, "\$x\$ (Non-Dim)")
ylabel!(p3, "\$|u|\$ (in m)")
title!(p3, "Displacement profiles for \$k_0 = "*string(k₀)*"\$ Nm\$^{-3}\$")

plt = plot(p1, p2, p3, layout=(3,1))
display(plt)
savefig(plt, "Example4.pdf")
