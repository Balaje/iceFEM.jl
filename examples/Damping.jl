using iceFEM
using Plots
using DelimitedFiles

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
plt2 = plot()

##################################
# Sample non-dimensionalization
##################################
ndproblem = non_dimensionalize(ice, fluid, 2π/100)

ωs = 2π*LinRange(0.01,0.125,400)

#########################################
# Obtain the Pierson Moskowitz spectrum
#########################################
ωₚ = 0.1*2π
ω₀ = 1.14*ωₚ
α = 0.0081
β₀ = 0.74
g = 9.8
S = (α*g^2 ./ωs.^5).*exp.(-β₀*(ω₀ ./ωs).^4) .+ 1e-10
strain²S = Vector{Float64}(undef,length(S))
fill!(strain²S,0.0)

strain_exp = readdlm("./examples/strain.txt",',', Float64)
plt3 = plot();
# plot!(plt3, strain_exp[:,1], strain_exp[:,2], label="")

α₀ = 0.1
#= A₀s = LinRange(1,40,400)
βs = reshape([(1 - A₀/(α₀ + 1im*ω)) for ω in ωs for A₀ in A₀s],(400,400))
pltC = contourf(ωs, A₀s, abs.(βs))
xlabel!(pltC, "\$ \\omega \$ (in Hz)")
ylabel!(pltC, "\$ |A_0| \$ (damping term)")
 =#

#for A₀ ∈ vcat(0,0.1,0.2,1,10:10:30)
for A₀ ∈ [0, 0.1, 0.2, 0.5, 1]
  Uₛ = zeros(Float64,length(ωs))
  ∂ₓ²Uₛ = zeros(Float64,length(ωs))
  for i in 1:lastindex(ωs)
    β = 1 - A₀/(α₀ + 1im*ωs[i])
    local fd = solve(ice, fluid, ωs[i], BeamType, WaterType; β = β)
    local x = 0:0.01:ndproblem.geo[1]
    Uₛ[i] = maximum(abs.(u₁(x, fd))) # Displacement
    ∂ₓ²Uₛ[i] = maximum(abs.(∂ₓ²u₁(x, fd)))*(ndproblem.geo[3])*(10^6) # Strain (in μ-strain)
    strain²S[i] = (∂ₓ²Uₛ[i])^2*S[i]
  end
  # Strain
  plot!(plt, ωs/(2π), ∂ₓ²Uₛ, label="Strain A₀ = "*string(A₀), line=(:solid, 2), yaxis=:log10)
  xlabel!(plt, "\$ \\omega/(2\\pi) \$ (in Hz)")
  ylabel!(plt, "\$ \\max \\,|\\partial_x^2 U| \\times 10^{-6} \$ ")
  # Displacement
  plot!(plt1, ωs/(2π), Uₛ, label="A₀ = "*string(A₀), line=(:solid, 2), yaxis=:log10)
  xlabel!(plt1, "\$ \\omega/(2\\pi) \$ (in Hz)")
  ylabel!(plt1, "\$ \\max \\,|U| \$ (in m)")
  # Strain Power spectral density
  if(10 ≤ A₀ ≤ 30)  
    plot!(plt2, ωs/(2π), strain²S, label="A₀ = "*string(A₀), line=(:solid, 3), yaxis=:identity)
  else
    plot!(plt2, ωs/(2π), strain²S, label="A₀ = "*string(A₀), line=(:dash, 2), yaxis=:identity)
  end 
  xlabel!(plt2, "\$ \\omega/(2\\pi) \$ (in Hz)")
  ylabel!(plt2, "\$ \\epsilon^2 \\,PSD \\times 10^{-12}\$")
  # xlabel!(plt3, "\$ \\omega/(2\\pi) \$ (in Hz)")
  # ylabel!(plt3, "Strain² \$\\times 10^{-12} \$")
end
plot!(plt2, strain_exp[:,1], strain_exp[:,2], label="Experimental", line=(:solid,2,:black), yaxis=:identity)

xlims!(plt, (0.01,0.125))
xlims!(plt1, (0.01,0.125))
xlims!(plt2, (0.01,0.125))
xlims!(plt3, (0.01,0.125))
#plt2 = plot(plt2, plt3, layout=(1,2))
plt = plot(plt, plt1, layout=(2,1))