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

mutable struct NonDimensionalProblem <: Any
  𝑙::Float64
  𝑘::Float64
  γ::Float64
  α::ComplexF64
  X::ComplexF64
  geo::Vector{Float64}
end

function preallocate_matrices(::Type{NonDimensionalProblem})
  NonDimensionalProblem(0, 0, 0, 0im, 0im, Vector{Float64}(undef,5))
end

function non_dimensionalize!(cache, Ice::Ice, Fluid::Fluid, ω)
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
  cache.𝑙 = (D/(ρₒ*g))^0.25
  cache.𝑘 = (k₀/(ρₒ*g))^0.25
  cache.γ = (ρᵢ/ρₒ)*(h/cache.𝑙)
  cache.α = ω^2*(cache.𝑙/g)

  d = cache.γ*cache.𝑙
  cache.X = (H-d)/(1im*ω*cache.𝑙^2)
  cache.geo[1] = L/cache.𝑙
  cache.geo[2] = H/cache.𝑙
  cache.geo[3] = h/cache.𝑙
  cache.geo[4] = x₀/cache.𝑙
  cache.geo[5] = g
  return nothing
end

function non_dimensionalize(Ice::Ice, Fluid::Fluid, ω)
  cache = preallocate_matrices(NonDimensionalProblem)
  non_dimensionalize!(cache, Ice, Fluid, ω)
  cache
end

## DataTypes describing the beam-type
struct FreeBedrock end
struct FreeClamped end
struct FreeHinged end
struct FreeFree end

## Function for the derivative of the displacements
