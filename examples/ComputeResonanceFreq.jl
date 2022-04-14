######################################################################
# Program to compute the resonant frequency of the fluid-ice system
######################################################################
using iceFEM
using Plots
using GenericLinearAlgebra
using LinearAlgebra

# Ice parameters
L = 3630
ρᵢ = 922.5
ν = 0.33
h = 280

ρₒ = 1025.0
g = 9.8
H = 500

function computeResonanceFrequency(ice, fluid, ω₀, verbosity=0)
  count = 1
  dw = 1e-6
  tol=1e-40
  fd = FiniteDepth(4)
  Δw = 1
  while (abs(Δw) > tol) && (count ≤ 50)
    s₁ = solve(ice, fluid, ω₀, FreeFree(), fd)
    s₂ = solve(ice, fluid, ω₀ + Δw, FreeFree(), fd)
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

Es = LinRange(1,5,30)
ω₀s = 2π*[0.022, 0.051]
ωₑs = [1/40, 1/15]
#ω₀s = 2π*[0.022, 0.051, big(0.12)]
#ωₑs = [1/40, 1/15, 1/10]
plt = plot()
clrs=[:red, :green]
ωᵣs = ones(Complex{Float64}, length(Es), length(ω₀s))
for (j,ω₀,ωₑ,clr) in zip(1:length(ω₀s), ω₀s,ωₑs,clrs)
  ωᵣs[:,j] = ones(Complex{Float64}, length(Es))*ω₀
  for i in 1:length(Es)
    ωᵣs[i,j] = computeResonanceFrequency(Ice(ρᵢ, Es[i]*1e9, ν, L, h),
                                         Fluid(ρₒ, 0, g, H, 0)
                                         ,(i==1) ? ω₀ : ωᵣs[i-1,j], 0)
  end
  plot!(plt, Es, real(ωᵣs[:,j])/(2π), linewidth=2, label="\$\\omega_r\$ (Theo)",
        linecolor=clr)
  plot!(plt, Es, ωₑ*Es.^0, linewidth=1, linestyle=:dash,
        linecolor=clr, label="\$T_{exp} = "*string(1/ωₑ)*"\$ s")
end
plot!(plt, Es, 0.12*Es.^0, linewidth=1, linestyle=:dash,
      linecolor=:blue, label="\$T_{exp} = "*string(round(1/0.12,digits=4))*"\$ s", legend=:topleft)
xlabel!(plt, "Young's Modulus \$E\$ (in GPa)")
ylabel!(plt, "\$\\omega_r/2\\pi\$ (in s\$^{-1}\$)")

plt1 = plot()
for (j, ω₀) in zip(1:length(ω₀s), ω₀s)
  plot!(plt1, Es, abs.(imag(ωᵣs[:,j]/(2π))), linewidth=2, label="Mode"*string(j), yaxis=:log)
end
xlabel!(plt1, "Young's Modulus \$E\$ (in GPa)")
ylabel!(plt1, "\$\\omega_i/2\\pi\$ (in s\$^{-1}\$)")


# Last one was using high precision arithmetic
omega_rs_mode3 = readdlm("omegars-mode3.txt", '\t', Complex{BigFloat}, '\n')
plot!(plt, Es, real(omega_rs_mode3/(2π)), linewidth=2, label="\$ \\omega_r\$ (Theo)",
      linecolor=:blue)
plot!(plt1, Es, abs.(imag(omega_rs_mode3/(2π))), linewidth=2, label="Mode"*string(3), yaxis=:log)
