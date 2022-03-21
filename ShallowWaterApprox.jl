include("Input.jl")

struct ShallowWater end
struct ShallowWaterSolution
  a₀::Vector{ComplexF64}
  b::Vector{ComplexF64}
  c::Vector{ComplexF64}
  p::Vector{ComplexF64}
  m::Vector{ComplexF64}
  ndp::NonDimensionalProblem
  BeamType
end
######################################################################################
# Function solve: Solves a frequency-domain problem
# Default option is to solve the ice-shelf problem using the FreeClamped() condition
#######################################################################################
function solve(Ice::Ice, Fluid::Fluid)
  solve_frequency_problem(Ice, Fluid, FreeClamped())
end

function solve(Ice::Ice, Fluid::Fluid, ω, ptype::Union{FreeClamped, FreeHinged}, ::ShallowWater)
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
  ShallowWaterSolution([a₀],0*c,c,0*m,m,ndp,ptype)
end
###########################################################
# u₁ is the displacement of the ice-shelf above the fluid #
###########################################################
function u₁(x, sol::ShallowWaterSolution)
  X = sol.ndp.X
  LL = sol.ndp.geo[1]
  m = sol.m
  c = sol.c
  @assert length(m) == length(c)
  if(sol.BeamType isa Union{FreeClamped, FreeHinged})
    return X*(c[1]*m[1]^2*exp.(m[1]*(x .- LL)) +
      c[2]*m[2]^2*exp.(m[2]*(x .- LL)) +
      c[3]*m[3]^2*exp.(m[3]*x) +
      c[4]*m[4]^2*exp.(m[4]*x) +
      c[5]*m[5]^2*exp.(m[5]*x) +
      c[6]*m[6]^2*exp.(m[6]*x))
  elseif(sol.BeamType isa FreeBedrock)
    X = sol.ndp.X
    xg = sol.ndp.geo[4]
    return X*(c[1]*m[1]^2*exp.(-m[1]*x) +
      c[2]*m[2]^2*exp.(-m[2]*x) +
      c[3]*m[3]^2*exp.(m[3]*x) +
      c[4]*m[4]^2*exp.(m[4]*x) +
      c[5]*m[5]^2*exp.(-m[5]*(x .- xg)) +
      c[6]*m[6]^2*exp.(-m[6]*(x .- xg)))
  end
end
#####################################
# Freq. Domain problem with bedrock #
#####################################
function solve(Ice::Ice, Fluid::Fluid, ω, ::FreeBedrock, ::ShallowWater)
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
  b = x[8:9]
  p = [p₁, p₂]
  ShallowWaterSolution([a₀], b, c, p, m, ndp, FreeBedrock())
end
#####################################################
# u₂ is the displacement profile above the bedrock
#####################################################
function u₂(x, sol::ShallowWaterSolution)
  xg = sol.ndp.geo[4]
  p = sol.p
  b = sol.b
  @assert length(p) == length(b)
  b[1]*exp.(-p[1]*(x .-xg)) + b[2]*exp.(-p[2]*(x .-xg))
end


## Solve the free-infinite problem
struct FreeInfinite end
function solve(Ice::Ice, Fluid::Fluid, ω, ::FreeInfinite, ::ShallowWater)
  ## To populate
  ω
end
function u₁(x, m, c, ndp::NonDimensionalProblem, ::FreeInfinite)
  ## To populate
  x
end

##############################
# Slope of the displacements #
##############################
function ∂ₓu₁(x, sol::ShallowWaterSolution)
  X = sol.ndp.X
  LL = sol.ndp.geo[1]
  m = sol.m
  c = sol.c
  𝑙 = sol.ndp.𝑙
  @assert length(m) == length(c)
  if(sol.BeamType isa Union{FreeClamped, FreeHinged})
    return (X/𝑙)*(c[1]*m[1]^3*exp.(m[1]*(x .- LL)) +
      c[2]*m[2]^3*exp.(m[2]*(x .- LL)) +
      c[3]*m[3]^3*exp.(m[3]*x) +
      c[4]*m[4]^3*exp.(m[4]*x) +
      c[5]*m[5]^3*exp.(m[5]*x) +
      c[6]*m[6]^3*exp.(m[6]*x))
  elseif(sol.BeamType isa FreeBedrock)
    X = sol.ndp.X
    xg = sol.ndp.geo[4]
    return (X/𝑙)*(-c[1]*m[1]^3*exp.(-m[1]*x) -
      c[2]*m[2]^3*exp.(-m[2]*x) +
      c[3]*m[3]^3*exp.(m[3]*x) +
      c[4]*m[4]^3*exp.(m[4]*x) -
      c[5]*m[5]^3*exp.(-m[5]*(x .- xg)) -
      c[6]*m[6]^3*exp.(-m[6]*(x .- xg)))
  end
end
function ∂ₓu₂(x, sol::ShallowWaterSolution)
  xg = sol.ndp.geo[4]
  p = sol.p
  b = sol.b
  @assert length(p) == length(b)
  (1/sol.ndp.𝑙)*(-p[1]*b[1]*exp.(-p[1]*(x .-xg))
                 -p[2]*b[2]*exp.(-p[2]*(x .-xg)))
end
