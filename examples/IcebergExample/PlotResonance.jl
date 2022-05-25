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
