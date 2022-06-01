using iceFEM
using Plots
using DelimitedFiles

ωᵣ = 2π*LinRange(0.01,0.125,300)
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

Es = LinRange(1,5,5)
L = 3630
ρᵢ = 922.5
Eᵢ = 2e9
ν = 0.33
h = 280

ρₒ = 1025.0
H = 500
fluid = Fluid(ρₒ, 0, g, H, 0)

ω₀s = [0.042, 0.056, 0.067, 0.076, 0.085]
ω₁s = ω₀s .+ 0.005

BeamType=FreeFree()
NModes=3
WaterType=FiniteDepth(NModes)

p = plot()
for (i,clr) in zip(1:length(Es), [:red,:green,:blue,:magenta,:black])
  local E = Es[i]
  ice = Ice(ρᵢ, E*1e9, ν, L, h)
  ∂ₓ²Uₛ = zeros(length(ωᵣ), 1)
  for j in 1:length(ωᵣ)
    fd = solve(ice, fluid, ωᵣ[j], BeamType, WaterType)
    x = 0:0.01:fd.ndp.geo[1]
    ∂ₓ²Uₛ[j] = maximum(abs.(∂ₓ²u₁(x,fd)))*(fd.ndp.γ*1/0.9)*(1/fd.ndp.𝑙)
  end
  strain²psd = ∂ₓ²Uₛ.^2 .*S
  plot!(p, ωᵣ/(2π), strain²psd, yaxis=:log10, label="\$ E = \$"*string(E)*" GPa",linecolor=clr)
  inds = (ωᵣ .> 2π*ω₀s[i]) .& (ωᵣ .< 2π*ω₁s[i])
  ωr1=ωᵣ[inds]/(2π)
  fωr1=strain²psd[inds]
  plot!(p, ωr1, fωr1, fillrange=1e-25*ωr1.^0,
        fillcolor=:red, fillalpha=0.5, label="")

  ### Find the area under the region and square root it
  ar = sum(fωr1[2:end-1]) + (fωr1[1]+fωr1[end])*0.5
  ar = ar*(ωᵣ[2]-ωᵣ[1]);
  ar = sqrt(ar)
  print("For E = "*string(E)*" GPa, RMS Strain = "*string(ar)*"\n")
end
#plot!(p,size=(800,500))
xlabel!(p,"\$ \\omega/(2\\pi) \$ (in \$s^{-1}\$)")
ylabel!(p,"\$ \\epsilon_{RMS} \$ ")
