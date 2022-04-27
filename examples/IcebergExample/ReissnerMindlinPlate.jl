using iceFEM
using Plots
using DelimitedFiles

ice = Ice(922.5, 5e9, 0.3, 3630, 280);
fluid = Fluid(1025, 0, 9.8, 500, 0);
μ = 10*(1+0.3)/(12+11*0.3)

##############################################
# Plot the roots of the RM-plate vs KL-plate #
##############################################
ωₛ = 2π*LinRange(0.01,0.125,200);
κₛₜ = zeros(ComplexF64, length(ωₛ), 3)
κₛᵣ = zeros(ComplexF64, length(ωₛ), 3)
for m=1:length(ωₛ)
  local ndp = non_dimensionalize(ice, fluid, ωₛ[m], ReissnerMindlinIce(); μ=μ)
  κₛₜ[m,:] = dispersion_ice(ndp.α, 1., ndp.γ, 2, ndp.geo[2]-ndp.γ)
  κₛᵣ[m,:] = dispersion_ice(ndp.α, 1., ndp.γ, ndp.geo[end], ndp.geo[3]^2, 2, ndp.geo[2]-ndp.γ)
end
ind = 1
plt = plot(ωₛ, abs.(κₛₜ[:,1]), label="KL-plate (Mode 1)", linewidth=2, linecolor="red")
plot!(plt, ωₛ, abs.(κₛᵣ[:,1]), label="RM-plate (Mode 1)", linewidth=2, linecolor="black")
plot!(plt, ωₛ, abs.(κₛₜ[:,2]), label="KL-plate (Mode 2)", linewidth=2, linecolor="red", linestyle=:dashdot)
plot!(plt, ωₛ, abs.(κₛᵣ[:,2]), label="RM-plate (Mode 2)", linewidth=2, linecolor="black", linestyle=:dashdot,legend=:topleft)

###########################################################
# Solve an example problem and verify energy conservation #
###########################################################
ω = 2π/15
sol = solve(ice, fluid, ω, FreeFree(), FiniteDepth(4), ReissnerMindlinIce(); μ=μ);
Aₚ = 9.8/(1im*ω)
R = sol.aₘ[1]/Aₚ; T = (sol.aₘ[6]/Aₚ)#*((κₜ*sinh(κₜ*(HH-γ)))/(kₜ*sinh(kₜ*(HH-γ))))
@show abs(R)^2 + abs(T)^2
x = LinRange(0,sol.ndp.geo[1],200);
Uᵣₘ = u₁(x, sol);
plt1 = plot(x, abs.(Uᵣₘ), label="R-M Plate", linewidth=2, color=:blue)
# Thin-plate solution
sol1 = solve(ice, fluid, ω, FreeFree(), FiniteDepth(4));
x = LinRange(0,sol1.ndp.geo[1],200);
Uₖₗ = u₁(x, sol1);
plot!(plt1, x, abs.(Uₖₗ), label="K-L Plate", linewidth=2, color=:red)
# Compare with LE
# ξ₂ = readdlm("./solDisp0_2GPa_10s_3630m_MS.dat", '\t', Float64, '\n')
# plot!(plt1, ξ₂[:,1], abs.(ξ₂[:,4]), color=:green, label="2D Elasticity",
#        linestyle=:dash,  linewidth=2)

########################################################
# Compute resonance frequency using the R-M Plate theory
########################################################
function computeResonanceFrequency(ice, fluid, ω₀, ::ReissnerMindlinIce, verbosity=0)
  count = 1
  dw = 1e-6
  tol=1e-10
  fd = FiniteDepth(4)
  Δw = 1
  while (abs(Δw) > tol) && (count ≤ 50)
    s₁ = solve(ice, fluid, ω₀, FreeFree(), fd, ReissnerMindlinIce())
    s₂ = solve(ice, fluid, ω₀ + Δw, FreeFree(), fd, ReissnerMindlinIce())
    Hⁿ = s₁.K
    Hⁿ⁺¹ = s₂.K

    ΔH = (Hⁿ⁺¹ - Hⁿ)/Δw
    dwₛ = eigvals(-ΔH\Hⁿ)
    (verbosity > 0) && print(string(Δw)*"\n")
    #dwₛ = eigvals(Hⁿ, -ΔH)
    dwₛ = sort(dwₛ,by=x->abs(x))
    Δw = dwₛ[1]
    #@time condH = cond(Hⁿ)
    (verbosity > 1) && print("Δw = "*string(Δw)*"\t ω₀ = "*string(ω₀)*"\t cond(H) = "*string(condH)*"...\n")
    if(!isinf(abs(Δw)))
      ω₀ = ω₀ + Δw
    end
    count+=1
  end
  (verbosity == 1) && print("Number of iterations = "*string(count)*"\n")
  ω₀
end

# Es = LinRange(1, 5, 30)
# ω₀s = 2π*[0.022, 0.051, 0.12]
# plt2 = plot()
# clrs = [:red, :green, :blue]
# ωᵣsᵣₘ = ones(ComplexF64, length(Es), length(ω₀s))
# ωᵣsₜₚ = ones(ComplexF64, length(Es), length(ω₀s))
# for (j,ω₀,clr) in zip(1:length(ω₀s), ω₀s, clrs)
#   for i in 1:length(Es)
#     ωᵣsₜₚ[i,j] = computeResonanceFrequency(Ice(922.5, Es[i]*1e9, 0.33, 3630, 280),
#                                            Fluid(1025, 0, 9.8, 500, 0)
#                                            ,(i==1) ? ω₀ : ωᵣsₜₚ[i-1,j], 0)
#     ωᵣsᵣₘ[i,j] = computeResonanceFrequency(Ice(922.5, Es[i]*1e9, 0.33, 3630, 280),
#                                            Fluid(1025, 0, 9.8, 500, 0)
#                                            ,(i==1) ? ω₀ : ωᵣsᵣₘ[i-1,j],
#                                            ReissnerMindlinIce(), 0)
#   end
#   plot!(plt2, Es, real(ωᵣsₜₚ[:,j])/(2π), linewidth=2, label="\$\\omega_r\$ (K-L)",
#         linecolor=clr, linestyle=:dash)
#   plot!(plt2, Es, real(ωᵣsᵣₘ[:,j])/(2π), linewidth=2, label="\$\\omega_r\$ (R-M)",
#         linecolor=clr, legend=:topleft)
# end
# xlabel!(plt2, "Young's Modulus \$E\$ (in GPa)")
# ylabel!(plt2, "\$\\omega_r/2\\pi\$ (in s\$^{-1}\$)")


# ωᵣLE1 = readdlm("./omegar_vs_youngs1_LE_real.txt", '\t', Float64, '\n')/(2π)
# ωᵣLE2 = readdlm("./omegar_vs_youngs2_LE_real.txt", '\t', Float64, '\n')/(2π)
# plot!(plt2, Es, ωᵣLE1, linewidth=1.5, label="\$\\omega_r\$ (LE)", linecolor=:red, linestyle=:dot)
# plot!(plt2, Es, ωᵣLE2, linewidth=1.5, label="\$\\omega_r\$ (LE)", linecolor=:green, linestyle=:dot)
