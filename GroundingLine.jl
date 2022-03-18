using Pkg
# Pkg.add("Polynomials")
using Polynomials
using Plots

# Physical quantities
ρₒ = 1025
ρᵢ = 922.5
Eᵢ = 2e9
ν = 0.3
k₀ = 1e8 # Spring constant for the bed.

# Dimensional quantities
L = 40000
H = 800
h = 500

# Derived quantities
d = (ρᵢ/ρₒ)*h
ω = 2π/2000
g = 9.8
D = Eᵢ*h^3/(12*(1-ν^2))
𝑙 = (D/(ρₒ*g))^0.25
𝑘 = (k₀/(ρₒ*g))^0.25
γ = d/𝑙
α = (ω^2*𝑙)/(g)
X = (H-d)/(1im*ω*𝑙^2)
Aₚ = g/(1im*ω)
xg = 0.7*L/𝑙# Grounding line position

# Solve the characteristic equation and collect the roots
# 1) Beam-Fluid
pl = Polynomial([α*𝑙/(H-d), 0, 1-γ*α, 0, 0, 0, 1])
m = roots(pl)
m = m[sortperm(real(m), rev=true)]

# 2) Beam-bedrock
#pl = Polynomial([1, 0, 0, 0, 1 + 𝑘^4 - γ*α])
pl = Polynomial([1 + 𝑘^4 - γ*α, 0, 0, 0, 1])
p = roots(pl)
p₁ = p[(real(p) .> 1e-9)][1]
p₂ = p[(real(p) .> 1e-9)][2]

# 3) Open-ocean
k = sqrt(α*𝑙/H)

# Construct the linear system
A = zeros(ComplexF64, 9, 9)
A[1,:] = [m[1]^4, m[2]^4, m[3]^4, m[4]^4, m[5]^4*exp(m[5]*xg), m[6]^4*exp(m[6]*xg), 0, 0, 0]
A[2,:] = [-m[1]^5, -m[2]^5, m[3]^5, m[4]^5, -m[5]^5*exp(m[5]*xg), -m[6]^5*exp(m[6]*xg), 0, 0, 0]
A[3,:] = [1, 1, 1, 1, exp(m[5]*xg), exp(m[6]*xg), -1, 0, 0]
A[4,:] = [-m[1], -m[2], m[3], m[4], -m[5]*exp(m[5]*xg), -m[6]*exp(m[6]*xg), 1im*k*H/(H-d), 0, 0]
A[5,:] = [X*m[1]^2*exp(-m[1]*xg), X*m[2]^2*exp(-m[2]*xg), X*m[3]^2*exp(m[3]*xg),
          X*m[4]^2*exp(m[4]*xg), X*m[5]^2, X*m[6]^2, 0,
          -1, -1]
A[6,:] = [-X*m[1]^3*exp(-m[1]*xg), -X*m[2]^3*exp(-m[2]*xg), X*m[3]^3*exp(m[3]*xg),
          X*m[4]^3*exp(m[4]*xg), -X*m[5]^3, -X*m[6]^3, 0,
          p₁, p₂]
A[7,:] = [X*m[1]^4*exp(-m[1]*xg), X*m[2]^4*exp(-m[2]*xg), X*m[3]^4*exp(m[3]*xg),
          X*m[4]^4*exp(m[4]*xg), X*m[5]^4, X*m[6]^4, 0,
          -p₁^2, -p₂^2]
A[8,:] = [-X*m[1]^5*exp(-m[1]*xg), -X*m[2]^5*exp(-m[2]*xg), X*m[3]^5*exp(m[3]*xg),
          X*m[4]^5*exp(m[4]*xg), -X*m[5]^5, -X*m[6]^5, 0,
          p₁^3, p₂^3]
A[9,:] = [-m[1]*exp(-m[1]*xg), -m[2]*exp(-m[2]*xg), m[3]*exp(m[3]*xg), m[4]*exp(m[4]*xg),
          -m[5], -m[6], 0, 0, 0]

f = [0, 0, Aₚ, 1im*k*Aₚ*H/(H-d), 0, 0, 0, 0, 0]

x = A\f

# Displacement profiles
c = x[1:6]
a₀ = x[7]
@show abs(a₀/Aₚ)
b = x[8:9]
p = [p₁, p₂]
function u₁(x, m, c)
  @assert length(m) == length(c)
  X*(c[1]*m[1]^2*exp.(-m[1]*x) +
    c[2]*m[2]^2*exp.(-m[2]*x) +
    c[3]*m[3]^2*exp.(m[3]*x) +
    c[4]*m[4]^2*exp.(m[4]*x) +
    c[5]*m[5]^2*exp.(-m[5]*(x .- xg)) +
    c[6]*m[6]^2*exp.(-m[6]*(x .- xg)))
end
function u₂(x, p, b)
  @assert length(p) == length(b)
  b[1]*exp.(-p[1]*(x .-xg)) + b[2]*exp.(-p[2]*(x .-xg))
end

x₁ = 0:0.01:xg

x₂ = xg:0.01:L/𝑙
U₁ = u₁(x₁, m, c)
U₂ = u₂(x₂, p, b)

plt = plot(x₁, abs.(U₁))
plot!(plt, x₂, abs.(U₂))
