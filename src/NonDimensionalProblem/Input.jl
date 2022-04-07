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

## DataTypes describing the beam-type
struct FreeBedrock end
struct FreeClamped end
struct FreeHinged end

## Function for the derivative of the displacements
