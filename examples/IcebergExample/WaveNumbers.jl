########################################
# Script to compare wave numbers
########################################
using iceFEM
using Plots

ice = Ice(922.5, 2e9, 0.3, 3630, 80)
fluid = Fluid(1025, 0, 9.8, 500, 0); # 500m water

########################################
# Plot the roots of the dispersion equation
########################################
ωₛ = 2π*LinRange(1/500,1/20,200);
NModes = 4
κₜ = zeros(ComplexF64, length(ωₛ),NModes+1)
κᵣ = zeros(ComplexF64, length(ωₛ),NModes+1)
k = zeros(ComplexF64, length(ωₛ),NModes+1)
R₁ = zeros(ComplexF64, length(ωₛ),NModes+1)
R₂ = zeros(ComplexF64, length(ωₛ),NModes+1)
R₃ = zeros(ComplexF64, length(ωₛ),NModes+1)
for m=1:length(ωₛ)
  local ndp = non_dimensionalize(ice, fluid, ωₛ[m], ReissnerMindlinIce())
  κₜ[m,:] = dispersion_ice(ndp.α, 1., ndp.γ, NModes, ndp.geo[2])/ndp.𝑙
  κᵣ[m,:] = dispersion_ice(ndp.α, 1., ndp.γ, ndp.geo[end]/0.9274, ndp.geo[3]^2, NModes, ndp.geo[2])/ndp.𝑙
  k[m,:] = dispersion_free_surface(ndp.α, NModes, ndp.geo[2])/ndp.𝑙

  R₁[m,:] = map(x -> iceFEM.DispersionEquations.f(x, ndp.α*ndp.geo[2],
                                                  1. /ndp.geo[2]^4,
                                                  ndp.γ/ndp.geo[2]),
                κₜ[m,:]*ndp.𝑙*1im*ndp.geo[2]) # Non-dim the input
  R₂[m,:] = map(x -> iceFEM.ReissnerMindlinPlate.f(x, ndp.α, 1., ndp.γ,
                                                   ndp.geo[end]/0.9274, ndp.geo[3]^2,
                                                   ndp.geo[2]),
                κₜ[m,:]*ndp.𝑙*1im)
  R₃[m,:] = map(x -> iceFEM.ReissnerMindlinPlate.f(x, ndp.α, 1., ndp.γ,
                                                   ndp.geo[end]/0.9274, ndp.geo[3]^2,
                                                   ndp.geo[2]),
                κᵣ[m,:]*ndp.𝑙*1im)
end
pl = plot(figsize=[20,4])
pl₂ = plot()
plₘ = plot()
#for ind = vcat(1:2,4:5)
for ind = [1,2,4,5]
  plot!(pl, 2π./ωₛ, ice.h^-1*2π./abs.(κₜ[:,ind]), linewidth=2, xlim=(20,200),
        yaxis=:log10, xaxis=:log10, label="(KL) Root = "*string(ind),
        linestyle=:dash)
  plot!(pl, 2π./ωₛ, ice.h^-1*2π./abs.(κᵣ[:,ind]), linewidth=2, xlim=(20,200),
        yaxis=:log10, xaxis=:log10, label="(RM) Root = "*string(ind),
        legend=:outertopright)
  xlabel!(pl, "\$ T \$ (in s)")
  ylabel!(pl, "\$ \\frac{\\lambda}{h} \$")

  ## Study difference in wave-lengths
  # plot!(pl₂, 2π./ωₛ, ice.L^-1*abs.(2π./abs.(κᵣ[:,ind]) - 2π./(abs.(κₜ[:,ind]))), linewidth=1,
  #       #yaxis=:log10, xaxis=:log10,
  #       label="H = "*string(fluid.H)*" m Mode = "*string(ind), legend=:topright)
  # xlabel!(pl₂, "\$ T \$ (in s)")
  # ylabel!(pl₂, "\$ \\frac{|\\lambda_{KL} - \\lambda_{RM}|}{L} \$")

  ## Study difference in wave-numbers
  plot!(plₘ, 2π./ωₛ, abs.(κᵣ[:,ind] - κₜ[:,ind])./abs.(κᵣ[:,ind]), linewidth=2,
        #yaxis=:log10, xaxis=:log10,
        #xlim=(20,200),
        label="Root = "*string(ind), legend=:topright)
  xlabel!(plₘ, "\$ T \$ (in s)")
  ylabel!(plₘ, "\$ |\\kappa_{KL} - \\kappa_{RM}|/|\\kappa_{RM}| \$")

  ## Study the residual
  plot!(pl₂, 2π./ωₛ, abs.(R₂[:,ind]), linewidth=2,
        #yaxis=:log10, xaxis=:log10,
        #xlim=(20,500),
        label="DEqᵣₘ(κₖₗ) Root = "*string(ind), legend=:topright)
  xlabel!(pl₂, "\$ T \$ (in s)")
  ylabel!(pl₂, "Residual")
end
title!(pl, "H = "*string(fluid.H)*" m")
title!(pl₂, "H = "*string(fluid.H)*" m")
title!(plₘ, "H = "*string(fluid.H)*" m")
vline!(pl, [50], label="T=50 s", linewidth=0.5, linestyle=:dash)
vline!(pl₂, [50], label="T=50 s", linewidth=0.5, linestyle=:dash)
vline!(plₘ, [50], label="T=50 s", linewidth=0.5, linestyle=:dash)

####### Solve the 50s problem
ω = 2π/12
sol = solve(ice, fluid, ω, FreeFree(), FiniteDepth(0), ReissnerMindlinIce(); μ=0.9274)
x = LinRange(0,sol.ndp.geo[1],200);
Uᵣₘ = u₁(x, sol);
pl₃ = plot(x, abs.(Uᵣₘ), label="R-M Plate", linewidth=2, color=:blue)
sol1 = solve(ice, fluid, ω, FreeFree(), FiniteDepth(0));
x = LinRange(0,sol1.ndp.geo[1],200);
Uₖₗ = u₁(x, sol1);
plot!(pl₃, x, abs.(Uₖₗ), label="K-L Plate", linewidth=2, color=:red)

######## Analyse the complex root with the thin-plate problem
fluid = Fluid(1025, 0, 9.8, 500, 0) # 500m water
hs = [12.5, 25, 50, 100, 200]
κᵣ = zeros(ComplexF64, length(ωₛ), NModes+1, length(hs))
κₜ = zeros(ComplexF64, length(ωₛ), NModes+1, length(hs))
for (h,j) in zip(hs,1:length(hs))
  local ice = Ice(922.5, 2e9, 0.3, 3630, h)
  for m in 1:length(ωₛ)
    local ndp = non_dimensionalize(ice, fluid, ωₛ[m], ReissnerMindlinIce())
    κₜ[m,:,j] = dispersion_ice(ndp.α, 1., ndp.γ, NModes, ndp.geo[2])/ndp.𝑙
    κᵣ[m,:,j] = dispersion_ice(ndp.α, 1., ndp.γ, ndp.geo[end]/0.9274, ndp.geo[3]^2, NModes, ndp.geo[2])/ndp.𝑙
  end
end

pls = fill(plot(),NModes+1)
for i in 1:NModes+1
  tmp = plot()
  for j in 1:length(hs)
    err = abs.(κᵣ[:,i,j]-κₜ[:,i,j])./abs.(κᵣ[:,i,j])
    (i == 5) ? plot!(tmp, 2π./ωₛ, err,  label="h = "*string(hs[j])*" m", linewidth=2, legend=:outertopright, xlim=(20,200)) :
      plot!(tmp, 2π./ωₛ, err, linewidth=2, label="",
            #xlim=(20,200)
            )
  end
  title!(tmp, "Root "*string(i))
  xlabel!(tmp, "\$ T \$ (in s)")
  ylabel!(tmp, "\$ |\\kappa_{KL} - \\kappa_{RM}|/|\\kappa_{RM}| \$")
  (i == 5) ? vline!(tmp, [50], label="T=50 s", linewidth=0.5, linestyle=:dash, legend=:topright) :
    vline!(tmp, [50], label="", linewidth=0.5, linestyle=:dash)
  pls[i] = tmp
end
pl₄ = plot(pls[1], pls[2], pls[4], pls[5], layout=(2,2))

########## Analyse the complex root with increasing μ
ice = Ice(922.5, 2e9, 0.3, 3630, 200)
fluid = Fluid(1025, 0, 9.8, 500, 0); # 500m water
ωₛ = 2π*LinRange(1/1000,1/400,200);
κₜ = zeros(ComplexF64, length(ωₛ),NModes+1)
κᵣ = zeros(ComplexF64, length(ωₛ),NModes+1)
κᵣ¹ = zeros(ComplexF64, length(ωₛ),NModes+1)
μ1 = 0.9274
pl₅ = plot()
for m=1:length(ωₛ)
  local ndp = non_dimensionalize(ice, fluid, ωₛ[m], ReissnerMindlinIce())
  κₜ[m,:] = dispersion_ice(ndp.α, 1., ndp.γ, NModes, ndp.geo[2])/ndp.𝑙
  κᵣ[m,:] = dispersion_ice(ndp.α, 1., ndp.γ, ndp.geo[end]/μ1, ndp.geo[3]^2, NModes, ndp.geo[2])/ndp.𝑙
end
μ2 = 8
for m=1:length(ωₛ)
  local ndp = non_dimensionalize(ice, fluid, ωₛ[m], ReissnerMindlinIce())
  κᵣ¹[m,:] = dispersion_ice(ndp.α, 1., ndp.γ, ndp.geo[end]/μ2, ndp.geo[3]^2, NModes, ndp.geo[2])/ndp.𝑙
end
for mode=[1,2,4]
  err = abs.(κᵣ[:,mode]-κₜ[:,mode])./abs.(κᵣ[:,mode])
  err1 = abs.(κᵣ¹[:,mode]-κₜ[:,mode])./abs.(κᵣ¹[:,mode])
  plot!(pl₅, 2π./ωₛ, err, label="|λₖ - λᵣ|/|λᵣ| Thick-plate μ = "*string(μ1)*" Mode "*string(mode), linewidth=1)
  plot!(pl₅, 2π./ωₛ, err1, label="|λₖ - λᵣ|/|λᵣ| Thick-plate μ = "*string(μ2)*" Mode "*string(mode), linewidth=1)#, xaxis=:log10, yaxis=:log10)
end
##############################
# Fox and Squire plot
##############################
# ice = Ice(922.5, 6e9, 0.3, 40000, 200)
# fluid = Fluid(1025, 0, 9.8, 1000, 0); #1000 m water (Replace 500m for the other depth)
# ωₛ = 2π*LinRange(1/60,1,200)
# κₜ = zeros(ComplexF64, length(ωₛ))
# κᵣ = zeros(ComplexF64, length(ωₛ))
# k = zeros(ComplexF64, length(ωₛ))
# #pl₁ = plot()
# for m=1:length(ωₛ)
#   local ndp = non_dimensionalize(ice, fluid, ωₛ[m], ReissnerMindlinIce())
#   κₜ[m] = dispersion_ice(ndp.α, 1., ndp.γ, 2, ndp.geo[2])[1]/ndp.𝑙
#   κᵣ[m] = dispersion_ice(ndp.α, 1., ndp.γ, ndp.geo[end]/0.9271, ndp.geo[3]^2, 2, ndp.geo[2])[1]/ndp.𝑙
#   k[m] = dispersion_free_surface(ndp.α, 2, ndp.geo[2])[1]/ndp.𝑙
# end
# plot!(pl₁, 2π./ωₛ, 2π./abs.(κₜ), linewidth=1, linecolor="red",
# yaxis=:log10, xaxis=:log10, label="1000 m (Thin-plate)", linestyle=:dot)

# plot!(pl₁, 2π./ωₛ, 2π./abs.(κᵣ), linewidth=1, linecolor="blue",
# yaxis=:log10, xaxis=:log10, label="1000 m (Thick-plate)", linestyle=:dot)

# plot!(pl₁, 2π./ωₛ, 2π./abs.(k), linewidth=1, linecolor="green",
# yaxis=:log10, xaxis=:log10, label="1000 m (Ocean)", linestyle=:dot,
# legend=:bottomright, xlim=(1,60))

#########################################################
# Plot the shear coefficient as a function of frequency
#########################################################
# ice = Ice(922.5, 2e9, 0.33, 40000, 200)
# fluid = Fluid(1025, 0, 9.8, 500, 0); # 100m water
# ωₛ = 2π*LinRange(0.01,0.1,200);
# μₛ = zeros(ComplexF64, length(ωₛ))
# function findμ(k, α, γ, δ, ζ, H)
#   den = (-α - k^5*tan(k*H) - γ*α*ζ/12*k^3*tan(k*H) - (1-γ*α)*k*tan(k*H))
#   num = (δ/ζ*(1-γ*α)*k^3*tan(k*H) + δ/12*(γ*α)^2*k*tan(k*H) - δ*α*γ/12*k*tan(k*H) - α*δ*k^2/ζ - γ*α^2*δ/12)
#   abs(num/den)
# end
# for m=1:length(ωₛ)
#   local ndp = non_dimensionalize(ice, fluid, ωₛ[m], ReissnerMindlinIce())
#   κₜ[m] = dispersion_ice(ndp.α, 1., ndp.γ, 2, ndp.geo[2]-ndp.γ)[1]
#   μₛ[m] = findμ(κₜ[m], ndp.α, ndp.γ, ndp.geo[end], ndp.geo[3]^2, ndp.geo[2] - ndp.γ)
# end
# pl₆ = plot(ωₛ/(2π), abs.(μₛ), linewidth=2, label="μ(ω)")
# xlabel!(pl₆, "\$\\omega/(2\\pi)\$ (Hz)")
# ylabel!(pl₆, "\$\\mu\$ (Shear correction factor)")
# for (clr,f) ∈ zip([:red, :green, :black], [1/20, 1/30, 1/50])
#   local ω = 2π*f
#   local ndp = non_dimensionalize(ice, fluid, ω, ReissnerMindlinIce())
#   k1 = dispersion_ice(ndp.α, 1., ndp.γ, 2, ndp.geo[2]-ndp.γ)[1]
#   μ1 = findμ(k1, ndp.α, ndp.γ, ndp.geo[end], ndp.geo[3]^2, ndp.geo[2] - ndp.γ)
#   @show μ1, f
#   vline!(pl₆, [f], label="T = "*string(1/f)*" s", linecolor=clr, linestyle=:dash, legend=false, linewidth=0.3)
#   hline!(pl₆, [μ1], label="μ = "*string(μ1), linecolor=clr, linestyle=:dash, legend=false, linewidth=0.3)
#   annotate!(pl₆, f + 0.2*f, μ1 + 0.2*μ1, text("(T,μ) = ("*string(1/f)*","*string(round(μ1, digits=3))*")", 10))
# end
