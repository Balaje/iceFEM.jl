using iceFEM
using Plots

##################################################
# Ice-shelf and fluid/bedrock parameters
##################################################
L = 40000
ρᵢ = 922.5
Eᵢ = 2e9
ν = 0.33
h = 500
ice = Ice(ρᵢ, Eᵢ, ν, L, h)

ρₒ = 1025.0
k₀ = 1e9
g = 9.8
H = 600
x₀ = 0.7*L
fluid = Fluid(ρₒ, k₀, g, H, x₀)

#################################################
# Eg 1: Solve one frequency domain problem
#################################################
BeamType = FreeBedrock()
WaterType = ShallowWater()
ω = 2π/200
sol₁ = solve(ice, fluid, ω, BeamType, WaterType)
Aₚ = g/(1im*ω)
@show abs(sol₁.a₀[1]/Aₚ)
# Visualize the solution
xg = sol₁.ndp.geo[4] # Grounding line location
LL = sol₁.ndp.geo[1] # Non dimensional length
𝑙 = abs(sol₁.ndp.𝑙)
x₁ = 0:0.01:xg; x₂ = xg:0.01:LL
U₁ = u₁(x₁, sol₁)
U₂ = u₂(x₂, sol₁)
plt1 = plot(x₁*𝑙, abs.(U₁), label="\$x < x_g\$", ylim=(0,3))
plot!(plt1, x₂*𝑙, abs.(U₂), label="\$ x > x_g\$", ylim=(0,3))
xlabel!(plt1, "\$x\$ (in m)")
ylabel!(plt1, "\$|u|\$ (in m)")
title!(plt1, "Displacement profiles for \$k_0 = "*string(k₀)*"\$ Nm\$^{-3}\$,\n Wave period = \$"*string(round(2π/ω, digits=4))*"\$ s")

#########################################################################
## Check displacement, slope, shear-force and bending-moment continuity #
#########################################################################
@assert u₁(xg, sol₁) ≈ u₂(xg, sol₁)
@assert ∂ₓu₁(xg, sol₁) ≈ ∂ₓu₂(xg, sol₁)
@assert ∂ₓ²u₁(xg, sol₁) ≈ ∂ₓ²u₂(xg, sol₁)
@assert ∂ₓ³u₁(xg, sol₁) ≈ ∂ₓ³u₂(xg, sol₁)

##########################################################################
## Check correctness of the non-local boundary condition (D2N map at GL) #
##########################################################################
𝐴 = (-sol₁.p[1])*(-sol₁.p[2])
𝐵 = (-sol₁.p[1] + -sol₁.p[2])
𝐶 = ( (-sol₁.p[1])^2 + (-sol₁.p[1])*(-sol₁.p[2]) + (-sol₁.p[2])^2 )
@assert -(𝐴)*(u₁(xg, sol₁)) + (𝐵)*(∂ₓu₁(xg,sol₁)) ≈ ∂ₓ²u₁(xg, sol₁)
@assert (-𝐴*𝐵)*(u₁(xg, sol₁)) + (𝐶)*(∂ₓu₁(xg,sol₁)) ≈ ∂ₓ³u₁(xg, sol₁)

##############################
# Eg 2: Displacement vs Freq.#
##############################
ωₛ = 2π*LinRange(0.001, 0.02, 500)
Uₛ¹ = zeros(length(ωₛ), 1)
Uₛ² = zeros(length(ωₛ), 1)

for i in 1:length(ωₛ)
  local sol = solve(ice, fluid, ωₛ[i], BeamType, WaterType)
  Uₛ¹[i] = maximum(abs.(u₁(x₁, sol)))
  Uₛ²[i] = maximum(abs.(u₂(x₂, sol)))
end
plt2 = plot(ωₛ, Uₛ¹, label="\$x < x_g\$")
plot!(plt2, ωₛ, Uₛ², label="\$x > x_g \$")
xlabel!(plt2,"\$\\omega\$ (in s\$^{-1}\$)")
ylabel!(plt2,"\$|u|\$ (in m)")
title!(plt2, "\$k_0 = "*string(k₀)*"\$ Nm\$^{-3}\$")

############################################################
# Eg 3: u(x₀) (Grounding Line disp.) vs k₀ (Spring Const.)
############################################################
k₀ₛ = 10 .^LinRange(6,9,100)
uₓ₀ₛ = zeros(length(k₀ₛ), 1)
∂ₓuₓ₀ₛ = zeros(length(k₀ₛ), 1)
plt3_1 = plot()
plt3_2 = plot()
for ω in [2π/50, 2π/100, 2π/200, 2π/8000]
  for i = 1:length(k₀ₛ)
    fl = Fluid(ρₒ, k₀ₛ[i], g, H, x₀)
    local sol = solve(ice, fl, ω, BeamType, WaterType)
    uₓ₀ₛ[i] = abs(u₁(sol.ndp.geo[4], sol)[1])
    ∂ₓuₓ₀ₛ[i] = abs(∂ₓu₁(sol.ndp.geo[4], sol)[1])
  end
  plot!(plt3_1, k₀ₛ, uₓ₀ₛ, label="\$T = \$ "*string(round((2π)/ω, digits=4))*" \$s\$")
  plot!(plt3_2, k₀ₛ, ∂ₓuₓ₀ₛ, label="\$T = \$ "*string(round((2π)/ω, digits=4))*" \$s\$")
end
xlabel!(plt3_1, "\$k_0\$ (Nm\$^{-3}\$)")
ylabel!(plt3_1, "\$|u(x_g)|\$")
ylabel!(plt3_2, "\$|\\partial_x u(x_g)|\$")
plt3 = plot(plt3_1,plt3_2,layout=(2,1))

############################################################
# Eg 3: u(x₀) (Grounding Line disp.) vs k₀ (Spring Const.)
############################################################
k₀ₛ = 10 .^LinRange(4,12,100)
uₓ₀ₛ = zeros(length(k₀ₛ), 1)
∂ₓuₓ₀ₛ = zeros(length(k₀ₛ), 1)
plt4_1 = plot()
plt4_2 = plot()
plt5_1 = plot()
plt5_2 = plot()

ω = 2π/200

𝐴s = zeros(ComplexF64, length(k₀ₛ), 1)
𝐵s = zeros(ComplexF64, length(k₀ₛ), 1)

for th in 200
  for i = 1:length(k₀ₛ)
    ic = Ice(ρᵢ, Eᵢ, ν, L, th)
    fl = Fluid(ρₒ, k₀ₛ[i], g, H, x₀)
    local sol = solve(ic, fl, ω, BeamType, WaterType)
    uₓ₀ₛ[i] = abs(u₁(sol.ndp.geo[4], sol)[1])
    ∂ₓuₓ₀ₛ[i] = abs(∂ₓu₁(sol.ndp.geo[4], sol)[1])

    r₁ = -sol.p[1]; r₂ = -sol.p[2];
    𝐴s[i] = abs(1/(r₁*r₂))
    𝐵s[i] = abs((r₁ + r₂)/(r₁*r₂))

  end
  plot!(plt4_1, k₀ₛ, uₓ₀ₛ, label="\$h = \$ "*string(round(th, digits=4))*" \$m\$")
  plot!(plt4_2, k₀ₛ, ∂ₓuₓ₀ₛ, label="\$h = \$ "*string(round(th, digits=4))*" \$m\$")

  plot!(plt5_1, k₀ₛ/10^6, abs.(𝐴s), xaxis=:log, legend = false, yguidefontsize=6)
  plot!(plt5_2, k₀ₛ/10^6, abs.(𝐵s), xaxis=:log, legend = false, yguidefontsize=6)
end
xlabel!(plt4_1, "\$k_0\$ (Nm\$^{-3}\$)")
ylabel!(plt4_1, "\$|u(x_g)|\$")
ylabel!(plt4_2, "\$|\\partial_x u(x_g)|\$")

xlabel!(plt5_1, "\$k_0\$ (MPa/m)")
xlabel!(plt5_2, "\$k_0\$ (MPa/m)")
ylabel!(plt5_1, "\$\\left|\\frac{1}{p_1p_2}\\right|\$")
ylabel!(plt5_2, "\$\\left|\\frac{p_1+p_2}{p_1p_2}\\right|\$")

plt4 = plot(plt4_1,plt4_2,layout=(2,1))
plt5 = plot(plt5_1,plt5_2,layout=(2,1))
title!(plt5, "Wave Period \$ T = "*string(round(2π/ω, digits=4))*"\$ s")
