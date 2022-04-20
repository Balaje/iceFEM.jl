using iceFEM
using Plots

#############################################
# Iceberg/fluid parameters
#############################################
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

ω = 2π/2000
#############################################
# Solve one frequency domain problem
#############################################
BeamType = FreeFree()
NModes = 3
WaterType = FiniteDepth(NModes)
ϕ = solve(ice, fluid, ω, BeamType, WaterType)
Aₚ = g/(1im*ω)
a₀ = ϕ.aₘ[1]/Aₚ
d₀ = ϕ.aₘ[NModes+2]/Aₚ
@show a₀ d₀
@show √(abs(a₀)^2 + abs(d₀)^2)

x = 0:0.01:ϕ.ndp.geo[1]
Uϕ = u₁(x, ϕ)
plt = plot(x,abs.(Uϕ), label="Displacement of Iceberg")
xlabel!(plt, "\$x\$ (Non-Dim)")
ylabel!(plt, "\$|u|\$ (in m)")

############################################
# Displacement of iceberg vs Frequency
###########################################
ωₛ = 2π*LinRange(0.01, 0.125, 400)
Uₛ = zeros(length(ωₛ), 1)
a₀ₛ = zeros(ComplexF64,length(ωₛ))
d₀ₛ = zeros(ComplexF64,length(ωₛ))
for i in 1:length(ωₛ)
  fd = solve(ice, fluid, ωₛ[i], BeamType, WaterType)
  Uₛ[i] = maximum(abs.(u₁(x,fd)))
  local Aₚ = g/(1im*ωₛ[i])
  a₀ₛ[i] = fd.aₘ[1]/Aₚ
  d₀ₛ[i] = fd.aₘ[NModes+2]/Aₚ
end
plt1 = plot(ωₛ/(2π), Uₛ, color=:red, linewidth=2)
xlabel!(plt1,"\$\\omega/2\\pi\$ (in s\$^{-1}\$)")
ylabel!(plt1,"\$|u|\$ (in m)")
plt2 = plot(ωₛ/(2π), abs.(a₀ₛ), color=:blue, linewidth=2, label="\$R(\\omega)\$")
plot!(plt2, ωₛ/(2π), abs.(d₀ₛ), color=:green, linewidth=2, label="\$T(\\omega)\$")
xlabel!(plt2,"\$\\omega/2\\pi\$ (in s\$^{-1}\$)")
ylabel!(plt2,"\$R(\\omega)\$ or \$T(\\omega)\$ (in m)")

savefig(plt1,"DispVsFreqEMM.pdf")
savefig(plt2,"RefTraVsFreqEMM.pdf")

################################################
# Reflection/Transmission vs Complex Frequency
################################################
# N = 200
# ωᵣ = 2π*LinRange(0.01,0.125,N)
# X = ωᵣ' .* ones(N)
# ωᵢ = LinRange(-0.08,0.08,N)
# Y = ones(N)' .* ωᵢ
# ωₛ = X + 1im*Y
# a₀ₛ = zeros(ComplexF64,N,N)
# d₀ₛ = zeros(ComplexF64,N,N)
# for i=1:size(ωₛ,1)
#   for j=1:size(ωₛ,2)
#     fd = solve(ice, fluid, ωₛ[i,j], BeamType, WaterType)
#     local Aₚ = g/(1im*ωₛ[i,j])
#     a₀ₛ[i,j] = fd.aₘ[1]/Aₚ
#     d₀ₛ[i,j] = fd.aₘ[NModes+2]/Aₚ
#   end
# end
# Rω = portrait(a₀ₛ,PTcgrid)
# Tω = portrait(d₀ₛ,PTcgrid)
# plt3 = plot(ωᵣ, ωᵢ, Rω)
