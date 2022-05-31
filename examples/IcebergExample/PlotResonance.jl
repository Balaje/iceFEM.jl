##################################################################
# Plot the resonance frequency as a function of Young's moduli
#         Run ResonanceVsYoungs.jl first
##################################################################
using Plots
using iceFEM
using DelimitedFiles

Es = LinRange(1,5,30)
ωₑs = [1/40, 1/15, 1/10]
ω₀s = 2π*[0.022, 0.051, big(0.12)]

clrs=[:red, :green, :black]

ωᵣs = readdlm("ResVsYoungs.txt", ',', Complex{BigFloat})
ωᵣs1 = readdlm("ResVsYoungs1.txt", ',', Complex{BigFloat})

plt = plot()
for (j, ω₀, ωₑ, clr) in zip(1:length(ω₀s), ω₀s, ωₑs, clrs)
  plot!(plt, Es, real(ωᵣs[:,j])/(2π), linewidth=2, label="\$\\omega_r\$ (L=3630 m)",
        linecolor=clr)
  plot!(plt, Es, ωₑ*Es.^0, linewidth=1, linestyle=:dashdot,
        linecolor=clr, label="\$T_{exp} = "*string(1/ωₑ)*"\$ s")
end
plot!(plt, Es, real(ωᵣs1[:,1])/(2π), linewidth=1, label="\$\\omega_r\$ (L=1740 m)",
      linecolor=:blue, linestyle=:dash)
plot!(plt, Es, 0.12*Es.^0, linewidth=1, linestyle=:dot,
      linecolor=:blue, label="\$T_{exp} = "*string(round(1/0.12,digits=4))*"\$ s",
      legend=:topleft)
xlabel!(plt, "Young's Modulus \$E\$ (in GPa)")
ylabel!(plt, "\$\\omega_r/2\\pi\$ (in s\$^{-1}\$)")

#0.013, 0.045
#0.05, 0.084
#0.088, 0.125
plot!(plt, Es, 0.013*Es.^0, fillrange=0.045*Es.^0, fillalpha=0.35, fillcolor=:red,label="")
plot!(plt, Es, 0.05*Es.^0, fillrange=0.084*Es.^0, fillalpha=0.35, fillcolor=:green, label="")
plot!(plt, Es, 0.088*Es.^0, fillrange=0.125*Es.^0, fillalpha=0.35, fillcolor=:blue, label="")


plt1 = plot()
for (j, ω₀) in zip(1:length(ω₀s), ω₀s)
  plot!(plt1, Es, abs.(imag(ωᵣs[:,j]/(2π))), linewidth=2,
        label="Mode"*string(j)*"(L=3630 m)", yaxis=:log)
end
plot!(plt1, Es, abs.(imag(ωᵣs1[:,1]/(2π))), linewidth=2,
      label="Mode"*string(1)*" (L=1740 m)", yaxis=:log, legend=:bottomleft)
xlabel!(plt1, "Young's Modulus \$E\$ (in GPa)")
ylabel!(plt1, "\$\\omega_i/2\\pi\$ (in s\$^{-1}\$)")

Y_err = sqrt.( abs.( (real(ωᵣs[:,1])/(2π) .- ωₑs[1]).^2
                     + (real(ωᵣs[:,2])/(2π) .- ωₑs[2]).^2
                     + (real(ωᵣs[:,3])/(2π) .- ωₑs[3]).^2 ) );
Y_err_1 = sqrt.( abs.(real(ωᵣs1[:,1])/(2π) .- 0.12).^2 );
plt3 = plot(Es, Y_err, linewidth=2, label="L = 3630 m")
plot!(plt3, Es, Y_err_1, linewidth=2, label="L = 1740 m")
xlabel!(plt3, "Young's Modulus \$E\$ (in GPa)")
ylabel!(plt3, "\$ \\sqrt{\\sum |\\omega_{r,theo} - \\omega_{exp}|^2}\$")

Y_err_2 = sqrt.(Y_err.^2 + Y_err_1.^2)
plot!(plt3, Es, Y_err_2, linewidth=2, linestyle=:dash,
      linecolor=:black, label="Overall Difference")
vline!(plt3, [Es[argmin(Y_err_2)]], linestyle=:dash,
       linewidth=1, linecolor=:green,
       label="\$E_{opt} = \$"*string(round(Es[argmin(Y_err_2)],digits=4))*" GPa")


########################################
# Find the Freq.vs Strain curve on Eopt
########################################
#L = 3630
Ls = [3630, 1740]
ρᵢ = 922.5
Eᵢ = 2e9
ν = 0.33
h = 280
g = 9.8
ρₒ = 1025.0
H = 500
Eopt = 1.9655
fluid = Fluid(ρₒ, 0, g, H, 0)
ωᵣ = 2π*vcat(LinRange(0.01,0.1125,500), LinRange(0.1125,0.125,2000))
plt4 = plot()
for (L,clr) in zip(Ls,[:red,:black])
  local ∂ₓ²Uₚ = zeros(length(ωᵣ), 1)
  local ice = Ice(ρᵢ, Eopt*1e9, ν, L, h)
  for j in 1:length(ωᵣ)
    fd = solve(ice, fluid, ωᵣ[j], FreeFree(), FiniteDepth(4))
    x = 0:0.01:fd.ndp.geo[1]
    local ∂ₓ²Uₚ[j] = maximum(abs.(∂ₓ²u₁(x,fd)))*(fd.ndp.γ*1/0.9)*(1/fd.ndp.𝑙)
  end
  #plt4 = plot(ωᵣ/(2π), ∂ₓ²Uₚ, linewidth=2, label="Theoretical strain E = 1.9655 GPa")
  plot!(plt4, ωᵣ/(2π), ∂ₓ²Uₚ, linewidth=2, label="E=1.9655 GPa, L="*string(L)*" m", linecolor=clr)
end
vspan!(plt4, [0.013,0.045], linecolor=:red, fillcolor=:red, fillalpha=0.5, label="")
vspan!(plt4, [0.05,0.084], linecolor=:green, fillcolor=:green, fillalpha=0.5,label="")
vspan!(plt4, [0.088,0.125], linecolor=:blue, fillcolor=:blue, fillalpha=0.5,label="",yaxis=:log10)
xlabel!(plt4, "\$\\omega/(2\\pi)\$ (in Hz)")
ylabel!(plt4, "\$ \\epsilon_{xx}\$")
