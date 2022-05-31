using iceFEM
using Plots

include("InterpolateFreq.jl");

ice = Ice(922.5, 1.9655e9, 0.3, 3630, 280)
fluid = Fluid(1025, 0, 9.8, 500, 0)

setprecision(BigFloat, 512)
ωᵒ = computeResonanceFrequency(ice, fluid, big(2π*0.12), 1);

Reωᵒ = real(ωᵒ)
Imωᵒ = imag(ωᵒ)

ωₚ = 0.092*2π
ω₀ = 1.14*ωₚ
α = 0.0081
β = 0.74
g = 9.8
S(ω) = (α*g^2 ./ω.^5).*exp.(-β*(ω₀ ./ω).^4)

# Compute damping
ϵₜₕ = √(S(Reωᵒ))
ϵₑₓ = 0.8e-6
Imωₑₓ = (ϵₜₕ/ϵₑₓ)^2*Imωᵒ

Q = 2π*0.1/Imωₑₓ
