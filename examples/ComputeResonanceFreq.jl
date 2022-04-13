######################################################################
# Program to compute the resonant frequency of the fluid-ice system
######################################################################
using iceFEM
using Plots
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
  tol=1e-12
  fd = FiniteDepth(4)
  Δw = 1
  while (abs(Δw) > tol) && (count ≤ 20)
    s₁ = solve(ice, fluid, ω₀, FreeFree(), fd)
    s₂ = solve(ice, fluid, ω₀ + Δw, FreeFree(), fd)
    Hⁿ = s₁.K
    Hⁿ⁺¹ = s₂.K

    ΔH = (Hⁿ⁺¹ - Hⁿ)/Δw
    dwₛ = eigvals(Hⁿ, -ΔH)
    dwₛ = sort(dwₛ,by=x->abs(x))
    Δw = dwₛ[1]
    condH = cond(Hⁿ)
    (verbosity > 0) && print("Δw = "*string(Δw)*"\t ω₀ = "*string(ω₀)*"\t det(H) = "*string(condH)*"...\n")
    if(!isinf(abs(Δw)))
      ω₀ = ω₀ + Δw
    end
    count+=1
  end
  ω₀
end


Es = LinRange(1,5,60)
ω₀s = 2π*[0.022, 0.051, 0.12]
ωₑs = [1/40, 1/15, 1/10]
plt = plot()
clrs=[:red, :green, :blue]
for (ω₀,ωₑ,clr) in zip(ω₀s,ωₑs,clrs)
  local ωᵣs = ones(ComplexF64, length(Es))*ω₀
  for i in 1:length(Es)
    local Eᵢ = Es[i]*1e9
    local ice = Ice(ρᵢ, Eᵢ, ν, L, h)
    local fluid = Fluid(ρₒ, 0, g, H, 0)
    ωᵣs[i] = computeResonanceFrequency(ice, fluid, (i==1) ? ω₀ : ωᵣs[i-1])
    #print("Real(ωᵣ/(2π)) = "*string((ωᵣs[i]/(2π)))*" For E = "*string(Eᵢ)*" GPa. \n")
  end
  plot!(plt, Es, real(ωᵣs)/(2π), linewidth=2, label="\$\\omega_r\$ (Theo)",
        linecolor=clr)
  plot!(plt, Es, ωₑ*Es.^0, linewidth=1, linestyle=:dash,
        linecolor=clr, label="\$T_{exp} = "*string(1/ωₑ)*"\$ s")
end
plot!(plt, Es, 0.12*Es.^0, linewidth=1, linestyle=:dash,
      linecolor=:blue, label="\$T_{exp} = "*string(round(1/0.12,digits=4))*"\$ s", legend=:outertopright)
xlabel!(plt, "Young's Modulus \$E\$ (in GPa)")
ylabel!(plt, "\$\\omega_r/2\\pi\$ (in s\$^{-1}\$)")
