using iceFEM
using Plots

###########################################
# Iceberg and fluid properties
###########################################
L = 3630
ρᵢ = 922.5
Eᵢ = 2.1034*10^9
ν = 0.33
h = 280
ice = Ice(ρᵢ, Eᵢ, ν, L, h)

ρ₀ = 1025
g = 9.8
H = 500
fluid = Fluid(ρ₀, 0, g, H, 0)

plot!(size=(400,200))

BeamType = FreeFree()
NModes = 3
WaterType = FiniteDepth(NModes)

plt = plot()
plt1 = plot()

##################################
# Sample non-dimensionalization 
##################################
ndproblem = non_dimensionalize(ice, fluid, 2π/100)

ωs = 2π*LinRange(0.01,0.125,400)

for A₀ ∈ [0, 0.001, 0.005, 0.025, 0.125]
Uₛ = zeros(length(ωs), 1)
∂ₓ²Uₛ = zeros(length(ωs), 1)
for i in 1:length(ωs)
    α₀ = A₀
    β = 1 - A₀/(α₀ + 1im*ωs[i])
    local fd = solve(ice, fluid, ωs[i], BeamType, WaterType; β = β)
    local x = 0:0.01:ndproblem.geo[1]
    Uₛ[i] = maximum(abs.(u₁(x, fd))) # Displacement 
    ∂ₓ²Uₛ[i] = maximum(abs.(∂ₓ²u₁(x, fd)))*(ndproblem.geo[3]) # Strain
end 

plot!(plt, ωs/2/π, ∂ₓ²Uₛ, label="Strain A₀ = "*string(A₀), yaxis=:log10, line=(:solid, 2))
xlabel!(plt, "\$ \\omega/(2\\pi) \$ (in Hz)")
ylabel!(plt, "\$ \\max \\,|\\partial_x^2 U| \$")

plot!(plt1, ωs/2/π, Uₛ, label="Displacement A₀ = "*string(A₀), yaxis=:log10, line=(:solid, 2))
xlabel!(plt1, "\$ \\omega/(2\\pi) \$ (in Hz)")
ylabel!(plt1, "\$ \\max \\,|U| \$ (in m)")
end