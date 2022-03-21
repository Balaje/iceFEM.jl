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
k₀ = 1e4
g = 9.8
H = 500
x₀ = 0.7*L
fluid = Fluid(ρₒ, k₀, g, H, x₀)

#################################################
# Eg 1: Solve one frequency domain problem
#################################################
BeamType = FreeBedrock()
WaterType = ShallowWater()
ω = 2π/4000
sol₁ = solve(ice, fluid, ω, BeamType, WaterType)
Aₚ = g/(1im*ω)
@show abs(sol₁.a₀[1]/Aₚ)
# Visualize the solution
xg = sol₁.ndp.geo[4] # Grounding line location
LL = sol₁.ndp.geo[1] # Non dimensional length
x₁ = 0:0.01:xg; x₂ = xg:0.01:LL
U₁ = u₁(x₁, sol₁)
U₂ = u₂(x₂, sol₁)
plt = plot(x₁, abs.(U₁), label="\$x < x_g\$")
plot!(plt, x₂, abs.(U₂), label="\$ x > x_g\$")
xlabel!(plt, "\$x\$ (Non-Dim)")
ylabel!(plt, "\$|u|\$ (in m)")
title!(plt, "Displacement profiles for \$k_0 = "*string(k₀)*"\$ Nm\$^{-3}\$")
display(plt); #readline()
savefig(plt, "Example1.pdf")

#################################################
# Eg 2: Displacement vs Freq.
#################################################
ωₛ = 2π*LinRange(0.001, 0.02, 500)
Uₛ¹ = zeros(length(ωₛ), 1)
Uₛ² = zeros(length(ωₛ), 1)
for i in 1:length(ωₛ)
  sol = solve(ice, fluid, ωₛ[i], BeamType, WaterType)
  Uₛ¹[i] = maximum(abs.(u₁(x₁, sol)))
  Uₛ²[i] = maximum(abs.(u₂(x₂, sol)))
end
plt = plot(ωₛ, Uₛ¹, label="\$x < x_g\$")
plot!(plt, ωₛ, Uₛ², label="\$x > x_g \$")
xlabel!(plt,"\$\\omega\$ (in s\$^{-1}\$)")
ylabel!(plt,"\$|u|\$ (in m)")
display(plt);
title!(plt, "\$k_0 = "*string(k₀)*"\$ Nm\$^{-3}\$")
savefig(plt, "Example2.pdf")

############################################################
# Eg 3: u(x₀) (Grounding Line disp.) vs k₀ (Spring Const.)
############################################################
k₀ₛ = 10 .^LinRange(6,9,100)
uₓ₀ₛ = zeros(length(k₀ₛ), 1)
∂ₓuₓ₀ₛ = zeros(length(k₀ₛ), 1)
p1 = plot()
p2 = plot()
for ω in [2π/50, 2π/100, 2π/200, 2π/8000]
  for i = 1:length(k₀ₛ)
    fl = Fluid(ρₒ, k₀ₛ[i], g, H, x₀)
    sol = solve(ice, fluid, ω, BeamType, WaterType)
    uₓ₀ₛ[i] = abs(u₁(sol.ndp.geo[4], sol))
    ∂ₓuₓ₀ₛ[i] = abs(∂ₓu₁(sol.ndp.geo[4], sol))
  end
  plot!(p1, k₀ₛ, uₓ₀ₛ, label="\$T = \$ "*string(round((2π)/ω, digits=4))*" \$s\$")
  plot!(p2, k₀ₛ, ∂ₓuₓ₀ₛ, label="\$T = \$ "*string(round((2π)/ω, digits=4))*" \$s\$")
end
xlabel!(p1, "\$k_0\$ (Nm\$^{-3}\$)")
ylabel!(p1, "\$|u(x_g)|\$")
ylabel!(p2, "\$|\\partial_x u(x_g)|\$")
plt = plot(p1,p2,layout=(2,1))
display(plt); #readline()
savefig(plt, "Example3.pdf")
