############################################################
# ... Continuing the wave-numbers script from WaveNumbers.jl
############################################################

using iceFEM
using Plots
ice = Ice(922.5, 2e9, 0.3, 3630, 80)
fluid = Fluid(1025, 0, 9.8, 500, 0); # 500m water
NModes = 4
κₜ = zeros(ComplexF64, length(ωₛ),NModes+1)
κᵣ = zeros(ComplexF64, length(ωₛ),NModes+1)
k = zeros(ComplexF64, length(ωₛ),NModes+1)
R₁ = zeros(ComplexF64, length(ωₛ),NModes+1)
R₂ = zeros(ComplexF64, length(ωₛ),NModes+1)
R₃ = zeros(ComplexF64, length(ωₛ),NModes+1)
μₛ = [1,2,4,8,16]
pl₇ = plot()
ωₛ = 2π*LinRange(1/500,1/20,200)
for μ in μₛ
  for m=1:length(ωₛ)
    local ndp = non_dimensionalize(ice, fluid, ωₛ[m], ReissnerMindlinIce())
    κₜ[m,:] = dispersion_ice(ndp.α, 1., ndp.γ, NModes, ndp.geo[2])/ndp.𝑙
    κᵣ[m,:] = dispersion_ice(ndp.α, 1., ndp.γ, ndp.geo[end]/μ, ndp.geo[3]^2, NModes, ndp.geo[2])/ndp.𝑙
    k[m,:] = dispersion_free_surface(ndp.α, NModes, ndp.geo[2])/ndp.𝑙

    R₁[m,:] = map(x -> iceFEM.DispersionEquations.f(x, ndp.α*ndp.geo[2],
                                                    1. /ndp.geo[2]^4,
                                                    ndp.γ/ndp.geo[2]),
                  κₜ[m,:]*ndp.𝑙*1im*ndp.geo[2]) # Non-dim the input
    R₂[m,:] = map(x -> iceFEM.ReissnerMindlinPlate.f(x, ndp.α, 1., ndp.γ,
                                                     ndp.geo[end]/μ, ndp.geo[3]^2,
                                                     ndp.geo[2]),
                  κₜ[m,:]*ndp.𝑙*1im)
    R₃[m,:] = map(x -> iceFEM.ReissnerMindlinPlate.f(x, ndp.α, 1., ndp.γ,
                                                     ndp.geo[end]/μ, ndp.geo[3]^2,
                                                     ndp.geo[2]),
                  κᵣ[m,:]*ndp.𝑙*1im)
  end

  choice = 2
  err = abs.(R₂[:,2]) # Only for the second root
  (choice==2) && (err = abs.(κᵣ[:,2] - κₜ[:,2])./abs.(κᵣ[:,2]))
  plot!(pl₇, 2π./ωₛ, err, legend=:topleft, label="μ = "*string(μ), linewidth=2)
  xlabel!(pl₇, "\$ T \$ (in s)")
  ylabel!(pl₇, "Residual")
  (choice==2) && ylabel!(pl₇, "Error")
end

#############################
# Plot μ as a function of ω #
#############################
function fμ(k, α, γ, δ, ζ, H)
  den = (α - k^5*tanh(k*H) + γ*α*ζ/12*k^3*tanh(k*H) - (1-γ*α)*k*tanh(k*H))
  num = (δ/ζ*(1-γ*α)*k^3*tanh(k*H) + δ/12*(γ*α)^2*k*tanh(k*H) - δ*α*γ/12*k*tanh(k*H) + α*δ*k^2/ζ + γ*α^2*δ/12)
  abs(num/den)
end
ice = Ice(922.5, 2e9, 0.3, 3630, 280)
pl₈ = plot()
ωₛ = 2π*LinRange(1/100,1/10,200)
μs = zeros(ComplexF64, length(ωₛ))
for m=1:length(ωₛ)
  local ndp = non_dimensionalize(ice, fluid, ωₛ[m], ReissnerMindlinIce())
  κₜ[m,:] = dispersion_ice(ndp.α, 1., ndp.γ, NModes, ndp.geo[2])
  μs[m] = fμ(1im*κₜ[m,1], ndp.α, ndp.γ, ndp.geo[end], ndp.geo[3]^2, ndp.geo[2])
end
plot!(pl₈, 2π./ωₛ, abs.(μs), label="", linewidth=2)
xlabel!(pl₈, "\$ T \$ (in s)")
ylabel!(pl₈, "\$ \\mu\$")

###################################################
# Plot the difference in the dispersion relations #
###################################################
function ΔDE(k, α, γ, δ, ζ, H)
  δ*((1-γ*α)*(γ*α/12 - k^2/ζ)*k*tan(k*H) + α*(k^2/ζ - γ)) + k*tan(k*H)*(γ*α*k^2*ζ)/12
end
ice = Ice(922.5, 2e9, 0.3, 3630, 200)
pl₉ = plot()
pl₁₀ = plot()
ωₛ = 2π*LinRange(1/1000,1/20,200)
ΔDₛ1 = zeros(ComplexF64, length(ωₛ))
ΔDₛ2 = zeros(ComplexF64, length(ωₛ))
αₛ = zeros(ComplexF64, length(ωₛ))
μ = 1
for m=1:length(ωₛ)
  local ndp = non_dimensionalize(ice, fluid, ωₛ[m], ReissnerMindlinIce())
  κₜ[m,:] = dispersion_ice(ndp.α, 1., ndp.γ, NModes, ndp.geo[2])
  κᵣ[m,:] = dispersion_ice(ndp.α, 1., ndp.γ, ndp.geo[end]/0.9274, ndp.geo[3]^2, NModes, ndp.geo[2])
  ΔDₛ1[m] = ΔDE(κₜ[m,1], ndp.α, ndp.γ, ndp.geo[end]/μ, ndp.geo[3]^2, ndp.geo[2])
  ΔDₛ2[m] = ΔDE(κₜ[m,2], ndp.α, ndp.γ, ndp.geo[end]/μ, ndp.geo[3]^2, ndp.geo[2])
  αₛ[m] = ndp.α/ndp.geo[2]
end
plot!(pl₉, 2π./ωₛ, abs.(ΔDₛ1), label="Travelling wave", linewidth=2)
plot!(pl₉, 2π./ωₛ, abs.(ΔDₛ2), label="Decaying wave", linewidth=2, xaxis=:log10, yaxis=:log10)
xlabel!(pl₉, "\$ T \$ (in s)")
ylabel!(pl₉, "\$ \\Delta DE\$")

plot!(pl₁₀, 2π./ωₛ, abs.(κₜ[:,1]), label="Travelling wave", linewidth=2)
plot!(pl₁₀, 2π./ωₛ, abs.(κᵣ[:,2]), label="Decaying wave", linewidth=2)
plot!(pl₁₀, 2π./ωₛ, abs.(sqrt.(αₛ)), label="|√(α(T)/H)|", linewidth=1, linestyle=:dash, linecolor=:black)
plot!(pl₁₀, 2π./ωₛ, ωₛ.^0, label="=1", linewidth=1, linestyle=:dashdotdot, linecolor=:black, xaxis=:log10, yaxis=:log10)
xlabel!(pl₁₀, "\$ T \$ (in s)")
ylabel!(pl₁₀, "\$ |\\kappa_{0,1}|\$")
