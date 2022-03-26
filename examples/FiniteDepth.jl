using iceFEM
using Plots

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
ω = 2π/500

#############################################################################
# Eg 1: Solve and compare finite depth solution  to shallow water solution  #
#############################################################################
BeamType = FreeClamped()
NModes = 3
fd = FiniteDepth(NModes)
solϕ = solve(ice, fluid, ω, BeamType, fd)
# Check the reflection coeffcients
Aₚ = g/(1im*ω)
@show (solϕ.aₘ[1]/Aₚ), abs(solϕ.aₘ[1]/Aₚ)
solSW = solve(ice, fluid, ω, BeamType, ShallowWater())
@show (solSW.a₀[1]/Aₚ), abs(solSW.a₀[1]/Aₚ)
# Check the profiles
x = 0:0.01:solϕ.ndp.geo[1]
Uϕ = u₁(x, solϕ)
USW = u₁(x, solSW)
plt = plot(x, abs.(Uϕ), label="Finite depth")
plot!(plt, x, abs.(USW), label="Shallow water")
xlabel!(plt, "\$x\$ (Non-Dim)")
ylabel!(plt, "\$|u|\$ (in m)")
title!(plt, "Clamped \$x = L\$, Wave Period = \$"*string(round(2π/ω, digits=4))*"\$ s")
display(plt)
