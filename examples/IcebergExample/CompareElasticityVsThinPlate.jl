###########################################################################
# Program to compare the elasticity solution with the thin-plate solution #
###########################################################################
using iceFEM
using DelimitedFiles
using Plots

# Load all the elasticity solution
L = 3630
ρᵢ = 922.5
Eᵢ = 2e9
ν = 0.33
h = 280
ice = Ice(ρᵢ, Eᵢ, ν, L, h)

ρₒ = 1025.0
g = 9.8
H = 500
fluid = Fluid(ρₒ, 0, g, H, 0)

ω = 2π/100
ϕ₁ = solve(ice, fluid, ω, FreeFree(), FiniteDepth(5))
x = 0:0.01:ϕ₁.ndp.geo[1]
Uϕ₁ = u₁(x, ϕ₁)
ξ₁ = readdlm("./solDisp0_2GPa_100s_3630m.dat", '\t', Float64, '\n')
plt = plot(x, abs.(Uϕ₁), color=:red, label="Thin Plate", linewidth=2)
plot!(plt, ξ₁[:,1], abs.(ξ₁[:,4]), color=:blue, label="2D Elasticity",
      linestyle=:dash,  linewidth=2)
xlabel!(plt, "\$\\bar{x}\$ (non Dim)")
ylabel!(plt, "\$|u|\$ (in m)")

ω = 2π/50
ϕ₂ = solve(ice, fluid, ω, FreeFree(), FiniteDepth(3))
x = 0:0.01:ϕ₂.ndp.geo[1]
Uϕ₂ = u₁(x, ϕ₂)
ξ₂ = readdlm("./solDisp0_2GPa_50s_3630m.dat", '\t', Float64, '\n')
plt1 = plot(x, abs.(Uϕ₂), color=:red, label="Thin Plate",  linewidth=2)
plot!(plt1, ξ₂[:,1], abs.(ξ₂[:,4]), color=:blue, label="2D Elasticity",
      linestyle=:dash,  linewidth=2)
xlabel!(plt1, "\$\\bar{x}\$ (non Dim)")
ylabel!(plt1, "\$|u|\$ (in m)")

############################################
# Displacement of iceberg vs Frequency
###########################################
#ωₛ = 2π*vcat(LinRange(0.01, 0.065,500), LinRange(0.065, 0.075, 1000),
#             LinRange(0.075, 0.117, 100), LinRange(0.118, 0.120, 4000))
#ωₛ = 2π*LinRange(0.119, 0.11925, 4000)
ωₛ = 2π*LinRange(0.01,0.125,1000)
Uₛ = zeros(length(ωₛ), 1)
a₀ₛ = zeros(ComplexF64,length(ωₛ))
d₀ₛ = zeros(ComplexF64,length(ωₛ))
for i in 1:length(ωₛ)
  fd = solve(ice, fluid, ωₛ[i], FreeFree(), FiniteDepth(3))
  local x = 0:0.01:fd.ndp.geo[1]
  Uₛ[i] = maximum(abs.(u₁(x,fd)))
end
plt2 = plot(ωₛ/(2π), Uₛ, color=:red, linewidth=2)
xlabel!(plt2,"\$\\omega/2\\pi\$ (in s\$^{-1}\$)")
ylabel!(plt2,"\$|u|\$ (in m)")
UₛLE = readdlm("./examples/dispVsFreq_2GPa_3630m.dat", '\t', Float64, '\n')
ωₛLE = 2π*LinRange(0.01,0.125,size(UₛLE,1))
# plot!(plt2, ωₛLE/(2π), UₛLE[:,2], color=:blue, linewidth=2, linestyle=:dash
#       ,xlim=[0.01,0.1], ylim=[0,2.5]
#       )
