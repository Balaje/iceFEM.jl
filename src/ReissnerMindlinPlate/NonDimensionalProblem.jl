function preallocate_matrices(::Type{NonDimensionalProblem}, ::Type{ReissnerMindlinIce})
  NonDimensionalProblem(0, 0, 0, 0im, 0im, Vector{Float64}(undef,6))
end

function non_dimensionalize!(cache, Ice::Ice, Fluid::Fluid, ω, ::Type{ReissnerMindlinIce}, μ)
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
  G = (Eᵢ/(2*(1+ν)))
  cache.geo[6] = (ρₒ*g*cache.γ)/(μ*G)
  return nothing
end

function non_dimensionalize(ice::Ice, fluid::Fluid, ω, ::ReissnerMindlinIce; μ=1)
  cache = preallocate_matrices(NonDimensionalProblem, ReissnerMindlinIce)
  non_dimensionalize!(cache, ice, fluid, ω, ReissnerMindlinIce, μ)
  cache
end
