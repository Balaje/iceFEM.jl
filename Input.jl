using Pkg
# Pkg.add("Polynomials")
using Polynomials
using Plots
using LaTeXStrings

struct Ice <: Any
  ρᵢ::Float64
  Eᵢ::Float64
  ν::Float64
  L::Float64
  h::Float64
end

struct Fluid <: Any
  ρₒ::Float64
  k₀::Float64
  g::Float64
  H::Float64
  x₀::Float64
end

struct NonDimensionalProblem <: Any
  𝑙::Float64
  𝑘::Float64
  γ::Float64
  α::ComplexF64
  X::ComplexF64
  geo::Vector{Float64}
end

function non_dimensionalize(Ice::Ice, Fluid::Fluid, ω)
  Eᵢ = Ice.Eᵢ
  ρᵢ = Ice.ρᵢ
  h = Ice.h
  ν = Ice.ν
  L = Ice.L

  ρₒ = Fluid.ρₒ
  k₀ = Fluid.k₀
  g = Fluid.g
  H = Fluid.H
  x₀ = Fluid.x₀

  D = Eᵢ*h^3/(12*(1-ν^2))
  𝑙 = (D/(ρₒ*g))^0.25
  𝑘 = (k₀/(ρₒ*g))^0.25
  γ = (ρᵢ/ρₒ)*(h/𝑙)
  α = ω^2*(𝑙/g)

  d = γ*𝑙
  X = (H-d)/(1im*ω*𝑙^2)
  geo = [L/𝑙, H/𝑙, h/𝑙, x₀/𝑙, g]
  NonDimensionalProblem(𝑙, 𝑘, γ, α, X, geo)
end

struct FreeBedrock end
function solve_frequency_problem(Ice::Ice, Fluid::Fluid, ω, ::FreeBedrock)
  ndp = non_dimensionalize(Ice, Fluid, ω)
  α = ndp.α
  𝑙 = ndp.𝑙
  g = ndp.geo[end]
  γ = ndp.γ
  𝑘 = ndp.𝑘
  X = ndp.X

  ω = √(α*g/𝑙)
  L = ndp.geo[1]*𝑙
  H = ndp.geo[2]*𝑙
  d = γ*𝑙
  Aₚ = g/(1im*ω)
  xg = ndp.geo[4]

  # Solve the beam-fluid dispersion equation
  pl = Polynomial([α*𝑙/(H-d), 0, 1-γ*α, 0, 0, 0, 1])
  m = roots(pl)
  m = m[sortperm(real(m), rev=true)]

  # Solve the beam-bedrock equation
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

  c = x[1:6]
  a₀ = x[7]
  #  @show abs(a₀/Aₚ)
  b = x[8:9]
  p = [p₁, p₂]

  c,a₀,b,m,p,ndp
end

function u₁(x, m, c, ndp::NonDimensionalProblem, ::FreeBedrock)
  X = ndp.X
  xg = ndp.geo[4]
  @assert length(m) == length(c)
  X*(c[1]*m[1]^2*exp.(-m[1]*x) +
    c[2]*m[2]^2*exp.(-m[2]*x) +
    c[3]*m[3]^2*exp.(m[3]*x) +
    c[4]*m[4]^2*exp.(m[4]*x) +
    c[5]*m[5]^2*exp.(-m[5]*(x .- xg)) +
    c[6]*m[6]^2*exp.(-m[6]*(x .- xg)))
end
function u₂(x, p, b, ndp::NonDimensionalProblem, ::FreeBedrock)
  @assert length(p) == length(b)
  xg = ndp.geo[4]
  b[1]*exp.(-p[1]*(x .-xg)) + b[2]*exp.(-p[2]*(x .-xg))
end

struct FreeClamped end
struct FreeHinged end
function solve_frequency_problem(Ice::Ice, Fluid::Fluid, ω, ptype::Union{FreeClamped, FreeHinged})
  ndp = non_dimensionalize(Ice, Fluid, ω)
  α = ndp.α
  𝑙 = ndp.𝑙
  g = ndp.geo[end]
  γ = ndp.γ
  𝑘 = ndp.𝑘
  X = ndp.X

  ω = √(α*g/𝑙)
  L = ndp.geo[1]*𝑙
  H = ndp.geo[2]*𝑙
  d = γ*𝑙
  Aₚ = g/(1im*ω)

  # 2) Ice-shelf
  pl = Polynomial([α*𝑙/(H-d), 0, 1-γ*α, 0, 0, 0, 1])
  m = roots(pl)
  m = m[sortperm(real(m), rev=true)]

  # 3) Open-ocean
  k = sqrt(α*𝑙/H)

  A = zeros(ComplexF64, 7, 7)
  A[1,:] = [m[1], m[2], m[3]*exp(m[3]*L/𝑙), m[4]*exp(m[4]*L/𝑙),
            m[5]*exp(m[5]*L/𝑙), m[6]*exp(m[6]*L/𝑙), 0]
  A[2,:] = [m[1]^2, m[2]^2, m[3]^2*exp(m[3]*L/𝑙), m[4]^2*exp(m[4]*L/𝑙),
            m[5]^2*exp(m[5]*L/𝑙), m[6]^2*exp(m[6]*L/𝑙), 0]
  if(ptype==FreeHinged())
    A[3,:] = [m[1]^4, m[2]^4, m[3]^4*exp(m[3]*L/𝑙), m[4]^4*exp(m[4]*L/𝑙),
              m[5]^4*exp(m[5]*L/𝑙), m[6]^4*exp(m[6]*L/𝑙), 0]
  else
    A[3,:] = [m[1]^3, m[2]^3, m[3]^3*exp(m[3]*L/𝑙), m[4]^3*exp(m[4]*L/𝑙),
              m[5]^3*exp(m[5]*L/𝑙), m[6]^3*exp(m[6]*L/𝑙), 0]
  end
  A[4,:] = [m[1]^4*exp(-m[1]*L/𝑙), m[2]^4*exp(-m[2]*L/𝑙),
            m[3]^4, m[4]^4, m[5]^4, m[6]^4, 0]
  A[5,:] = [m[1]^5*exp(-m[1]*L/𝑙), m[2]^5*exp(-m[2]*L/𝑙),
            m[3]^5, m[4]^5, m[5]^5, m[6]^5, 0]
  A[6,:] = [exp(-m[1]*L/𝑙), exp(-m[2]*L/𝑙), 1, 1, 1, 1, -1]
  A[7,:] = [m[1]*exp(-m[1]*L/𝑙), m[2]*exp(-m[2]*L/𝑙), m[3],
            m[4], m[5], m[6], 1im*k*H/(H-d)]

  b = [0, 0, 0, 0, 0, Aₚ, Aₚ*1im*k*H/(H-d)]

  x = A\b
  c = x[1:6]
  a₀ = x[7]
  c,a₀,m,ndp
end
function u₁(x, m, c, ndp::NonDimensionalProblem, ::Union{FreeClamped, FreeHinged})
  X = ndp.X
  LL = ndp.geo[1]
  X*(c[1]*m[1]^2*exp.(m[1]*(x .- LL)) +
    c[2]*m[2]^2*exp.(m[2]*(x .- LL)) +
    c[3]*m[3]^2*exp.(m[3]*x) +
    c[4]*m[4]^2*exp.(m[4]*x) +
    c[5]*m[5]^2*exp.(m[5]*x) +
    c[6]*m[6]^2*exp.(m[6]*x)
     )
end

## Function for the derivative of the displacements
function ∂ₓu₁(x, m, c, ndp::NonDimensionalProblem, ::FreeBedrock)
  X = ndp.X
  xg = ndp.geo[4]
  𝑙=ndp.𝑙
  @assert length(m) == length(c)
  (X/𝑙)*(-c[1]*m[1]^3*exp.(-m[1]*x) -
    c[2]*m[2]^3*exp.(-m[2]*x) +
    c[3]*m[3]^3*exp.(m[3]*x) +
    c[4]*m[4]^3*exp.(m[4]*x) -
    c[5]*m[5]^3*exp.(-m[5]*(x .- xg)) -
    c[6]*m[6]^3*exp.(-m[6]*(x .- xg)))
end
function ∂ₓu₁(x, m, c, ndp::NonDimensionalProblem, ::Union{FreeClamped, FreeHinged})
  X = ndp.X
  LL = ndp.geo[1]
  (X/LL)*(c[1]*m[1]^3*exp.(m[1]*(x .- LL)) +
    c[2]*m[2]^3*exp.(m[2]*(x .- LL)) +
    c[3]*m[3]^3*exp.(m[3]*x) +
    c[4]*m[4]^3*exp.(m[4]*x) +
    c[5]*m[5]^3*exp.(m[5]*x) +
    c[6]*m[6]^3*exp.(m[6]*x)
     )
end

struct FreeInfinite end
function solve_frequency_problem(Ice::Ice, Fluid::Fluid, ω, ::FreeInfinite)
  ## To populate
  ω
end
function u₁(x, m, c, ndp::NonDimensionalProblem, ::FreeInfinite)
  ## To populate
  x
end

## Solve the free-infinite problem
