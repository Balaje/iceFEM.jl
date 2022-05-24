using iceFEM
using Plots
using ComplexPhasePortrait

#############################################
# Iceberg/fluid parameters
#############################################
L = 3630
ρᵢ = 922.5
Eᵢ = 2e9
ν = 0.33
h = 280
ice = Ice(ρᵢ, Eᵢ, ν, L, h)

ρₒ = 1025.0
g = 9.8
H = 500
fluid = Fluid(ρₒ, 0, g, H, 0)

##########################################################################
# Solve a set of frequency domain problems in the real frequency space
##########################################################################
BeamType=FreeFree()
NModes=3
WaterType=FiniteDepth(NModes)
# ωᵣ = 2π*LinRange(0.01,0.125,200)
# Ω₁ = FreqSpace{Float64}(4*NModes+6, length(ωᵣ)) # Real Coarse-Frequency space
# for i in 1:length(ωᵣ)
#   fd = solve(ice, fluid, ωᵣ[i], BeamType, WaterType)
#   Ω₁.ω[i] = ωᵣ[i]
#   Ω₁.H₀[i] = fd.K
#   Ω₁.f₀[i] = reshape(fd.f,length(fd.f),1)
# end
# ω₁s,H₁s,f₁s = InterpolateFreqDomain(Ω₁, (1600,));
# for i in 1:length(f₁s)
#   Aₚ = 9.8/(1im*ω₁s[i])
#   λ = (H₁s[i])\(f₁s[i])
#   a₀ = λ[1]/Aₚ
#   d₀ = λ[3NModes+8]/Aₚ
# end

# #############
# # Visualize #
# #############
# plt = plot(ωᵣ, [abs.((Ω₁.H₀[i]\Ω₁.f₀[i]))[1]*ωᵣ[i]/9.8 for i in 1:length(ωᵣ)],
#            linecolor=:red, linestyle=:dash, label="Coarse Solution")
# plot!(plt, ω₁s, [abs.((H₁s[i]\f₁s[i]))[1]*ω₁s[i]/9.8 for i in 1:length(ω₁s)],
#      linecolor=:blue, linestyle=:solid, label="Refined Solution")

###########################################################################
# Solve a set of frequency domain problems in the complex frequency space #
###########################################################################
xc = 2π*LinRange(0.01,0.04,50)
yc = LinRange(-0.06,0.06,50)
ωc = (xc' .* ones(length(yc))) + 1im*(ones(length(xc))' .* yc)
# Ω₂ = FreqSpace{ComplexF64}(4*NModes+6, size(ωc,1))
# for i in 1:size(ωc,1)
#   for j in 1:size(ωc,2)
#     fd = solve(ice, fluid, ωc[i,j], BeamType, WaterType)
#     Ω₂.ω[i,j] = ωc[i,j]
#     Ω₂.H₀[i,j] = fd.K
#     Ω₂.f₀[i,j] = reshape(fd.f,length(fd.f),1)
#   end
# end
# ω₂s,H₂s,f₂s = InterpolateFreqDomain(Ω₂, (300,300));
#############
# Visualize #
#############
# a₀s = zeros(ComplexF64,size(ω₂s))
# d₀s = zeros(ComplexF64,size(ω₂s))
# for i in 1:size(a₀s,1)
#   for j in 1:size(a₀s,2)
#     Aₚ = 9.8/(1im*ω₂s[i,j])
#     a₀s[i,j] = ((H₂s[i,j])\(f₂s[i,j]))[1]/Aₚ
#     d₀s[i,j] = ((H₂s[i,j])\(f₂s[i,j]))[3NModes+8]/Aₚ
#   end
# end
Rω = portrait(a₀s, PTproper, ctype="nist");
Tω = portrait(d₀s, PTproper, ctype="nist");
plt1 = plot(LinRange(0.01,0.04,300), LinRange(-0.06,0.06,300)/(2π), reverse(Rω,dims=1))
xlims!(plt1, (0.01,0.04))
ylims!(plt1, (-0.06,0.06)./(2π))
plt2 = plot(LinRange(0.01,0.04,300), LinRange(-0.06,0.06,300)/(2π), reverse(Tω,dims=1))
xlims!(plt2, (0.01,0.04))
ylims!(plt2, (-0.06,0.06)./(2π))

ω₀s = zeros(ComplexF64,2,1)
ω₀s[1] = computeResonanceFrequency(ice, fluid, 2π*0.012, 1)
ω₀s[2] = computeResonanceFrequency(ice, fluid, 2π*0.022, 1)
scatter!(plt1, real(ω₀s)/(2π), imag(ω₀s)/(2π), markercolor=[:white,:white], legend=false)
scatter!(plt2, real(ω₀s)/(2π), imag(ω₀s)/(2π), markercolor=[:white,:white], legend=false)

xlabel!(plt1, "\$ Re(\\omega)/(2\\pi)\$ (in \$s^{-1}\$)")
ylabel!(plt1, "\$ Im(\\omega)/(2\\pi)\$ (in \$s^{-1}\$)")
xlabel!(plt2, "\$ Re(\\omega)/(2\\pi)\$ (in \$s^{-1}\$)")
ylabel!(plt2, "\$ Im(\\omega)/(2\\pi)\$ (in \$s^{-1}\$)")
