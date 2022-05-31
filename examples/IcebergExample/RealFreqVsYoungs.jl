using iceFEM
using Plots
using DelimitedFiles
using FindPeaks1D
using Images

#ωᵣ = 2π*LinRange(0.01,0.125,300)
ωᵣ = 2π*LinRange(0.0375,0.1,300)
#ωᵣ = 2π*LinRange(0.01,0.0375,300)
######################################################
# First obtain the Pierson Moskowitz Spectrum
######################################################
ωₚ = 0.092*2π
ω₀ = 1.14*ωₚ
α = 0.0081
β = 0.74
g = 9.8
S=(α*g^2 ./ωᵣ.^5).*exp.(-β*(ω₀ ./ωᵣ).^4) .+ 1e-15
S_exp = readdlm("inc_wave.txt",',', Float64)

Es = LinRange(1,5,5) # And do it with 600 points
ε̄ = zeros(Float64, length(Es))
L = 3630
ρᵢ = 922.5
Eᵢ = 2e9
ν = 0.33
h = 280

ρₒ = 1025.0
H = 500
fluid = Fluid(ρₒ, 0, g, H, 0)

BeamType=FreeFree()
NModes=3
WaterType=FiniteDepth(NModes)

p = plot()

clrs = repeat([:red,:green,:blue,:magenta,:black], length(Es))
for (i,clr) in zip(1:length(Es), clrs)
  local E = Es[i]
  ice = Ice(ρᵢ, E*1e9, ν, L, h)
  ∂ₓ²Uₛ = zeros(length(ωᵣ), 1)
  for j in 1:length(ωᵣ)
    fd = solve(ice, fluid, ωᵣ[j], BeamType, WaterType)
    x = 0:0.01:fd.ndp.geo[1]
    ∂ₓ²Uₛ[j] = maximum(abs.(∂ₓ²u₁(x,fd)))*(fd.ndp.γ*1/0.9)*(1/fd.ndp.𝑙)
  end
  strain²psd = ∂ₓ²Uₛ.^2 .*S
  pkindices,  = findpeaks1d(vec(strain²psd), prominence=1e-17, width=(0,5))
  # pkindices,  = findpeaks1d(vec(strain²psd), prominence=1e-22)
  scatter!(p, ωᵣ[pkindices]/(2π), strain²psd[pkindices],
           color="red", markersize=5, label="")

  ######################################
  # Highlight the area under the peaks #
  ######################################
  δ⁻w = pkindices .- 15
  δ⁺w = pkindices .+ 15
  inds = vcat(δ⁻w[1], δ⁺w[1])
  ωᵣ1 = ωᵣ[inds[1]:inds[2]]/(2π)
  fωᵣ1 = strain²psd[inds[1]:inds[2]]

  plot!(p, ωᵣ1, fωᵣ1, fillrange=1e-25*ωᵣ1.^0,
        fillcolor=:red, fillalpha=0.5, label="", legend=:bottomright)
  plot!(p, ωᵣ/(2π), strain²psd,
        yaxis=:log10, label="\$ E = \$"*string(E)*" GPa",
        linecolor=clr)

  #######################################################
  # Find the area under the region and square root it   #
  #######################################################
  ε̄[i] = sum(fωᵣ1[2:end-1]) + (fωᵣ1[1]+fωᵣ1[end])*0.5
  ε̄[i] = ε̄[i]*(ωᵣ[2]-ωᵣ[1]);
  ε̄[i] = sqrt(ε̄[i])
  #print("For E = "*string(E)*" GPa, RMS Strain = "*string(ε̄[i])*"\n")
end
xlabel!(p,"\$ \\omega/(2\\pi) \$ (in \$s^{-1}\$)")
ylabel!(p,"\$ \\epsilon^2 PSD \$ ")

######################
# Compare RMS Strain #
######################
pp = plot(Es, ε̄, label="RMS Strain")
plot!(pp, Es, 4e-6*Es.^0, label="\$\\varepsilon_{RMS,Exp} = 4\\,\\mu \$ strain",yaxis=:log10)
