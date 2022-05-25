######################################################################
# Program to compute the resonant frequency of the fluid-ice system
######################################################################
using iceFEM
using Plots
using GenericLinearAlgebra
using LinearAlgebra


# Arbitrary precision arithmetic
setprecision(BigFloat, 512)

# Ice parameters
L = 1740
ρᵢ = 922.5
ν = 0.33
h = 280
ρₒ = 1025.0
g = 9.8
H = 500

Es = LinRange(1,5,30)
#ω₀s = 2π*[0.022, 0.051, big(0.12)]
#ωₑs = [1/40, 1/15, 1/10]
ω₀s = 2π*[big(0.12)]
ωₑs = [0.12]
plt = plot()
clrs=[:red, :green, :black]
ωᵣs = ones(Complex{BigFloat}, length(Es), length(ω₀s))
for (j,ω₀,ωₑ,clr) in zip(1:length(ω₀s), ω₀s,ωₑs,clrs)
  ωᵣs[:,j] = ones(Complex{BigFloat}, length(Es))*ω₀
  for i in 1:length(Es)
    ωᵣs[i,j] = computeResonanceFrequency(Ice(ρᵢ, Es[i]*1e9, ν, L, h),
                                         Fluid(ρₒ, 0, g, H, 0)
                                         ,(i==1) ? ω₀ : ωᵣs[i-1,j], 1)
  end
end

########################################################################
# Command to read and write files containing the resonance frequencies #
########################################################################
# open("ResVsYoungs.txt", "w") do io
#   writedlm(io, ωᵣs, ',')
# end;
# readdlm("ResVsYoungs.txt", ',', Complex{BigFloat})
