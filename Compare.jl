include("Input.jl")

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

ω = 2π/1000
#################################################
# Eg 1: Solve one frequency domain problem for clamped ice-shelf
#################################################
c,a₀,m,ndp = solve_frequency_problem(ice, fluid, ω, FreeClamped())
Aₚ = ndp.geo[end]/(1im*ω)
@show abs(a₀/Aₚ)
# Visualize the solution
x₁ = 0:0.05:ndp.geo[1]
U₁ = u₁(x₁, m, c, ndp, FreeClamped())
p1 = plot(x₁, abs.(U₁), legend=false)
xlabel!(p1, "\$x\$ (Non-Dim)")
ylabel!(p1, "\$|u|\$ (in m)")
title!(p1, "Displacement profile")

#################################################
# Eg 1: Solve one frequency domain problem for hinged ice-shelf
#################################################
c,a₀,m,ndp = solve_frequency_problem(ice, fluid, ω, FreeHinged())
Aₚ = ndp.geo[end]/(1im*ω)
@show abs(a₀/Aₚ)
# Visualize the solution
x₁ = 0:0.05:ndp.geo[1]
U₁ = u₁(x₁, m, c, ndp, FreeClamped())
p2 = plot(x₁, abs.(U₁), legend=false)
xlabel!(p2, "\$x\$ (Non-Dim)")
ylabel!(p2, "\$|u|\$ (in m)")
title!(p2, "Displacement profile")

#################################################
# Eg 3: Solve one frequency domain problem for bedrock ice-shelf
# Note: Now the length should change along with the bedrock parameters
#################################################
L = 40000/0.7
ice = Ice(ρᵢ, Eᵢ, ν, L, h)
k₀ = 1e6 # Stiffness of the bedrock
x₀ = 0.7*L # For the grounding line problem
fluid = Fluid(ρₒ, k₀, g, H, x₀)
c,a₀,b,m,p,ndp = solve_frequency_problem(ice, fluid, ω, FreeBedrock())
Aₚ = ndp.geo[end]/(1im*ω)
@show abs(a₀/Aₚ)
# Visualize the solution
x₁ = 0:0.01:ndp.geo[4]
U₁ = u₁(x₁, m, c, ndp, FreeBedrock())
p3 = plot(x₁, abs.(U₁), label="\$x < x_g\$")
x₂ = ndp.geo[4]:0.01:ndp.geo[1]
U₂ = u₂(x₂, p, b, ndp, FreeBedrock())
plot!(p3, x₂, abs.(U₂), label="\$ x > x_g\$")
xlabel!(p3, "\$x\$ (Non-Dim)")
ylabel!(p3, "\$|u|\$ (in m)")
title!(p3, "Displacement profiles for \$k_0 = "*string(k₀)*"\$ Nm\$^{-3}\$")

plt = plot(p1, p2, p3, layout=(3,1))
display(plt)
savefig(plt, "Example.pdf")
