struct ShallowWater end
mutable struct ShallowWaterSolution{T<:ComplexF64}
  a₀::Vector{T}
  b::Vector{T}
  c::Vector{T}
  p::Vector{T}
  m::Vector{T}
  ndp::NonDimensionalProblem
  A::Matrix{T}
  x::Vector{T}
  BeamType
end

#### High level-interface
function solve(Ice::Ice, Fluid::Fluid, ω, ptype, ::ShallowWater)
  cache = preallocate_matrices(ShallowWater(), ptype)
  solve!(cache, Ice, Fluid, ω, ptype, ShallowWater())
  cache[1]
end
function u₁(x, sol::ShallowWaterSolution)
  op = Vector{ComplexF64}(undef, length(x))
  u₁!(op, x, sol)
  op
end
function ∂ₓu₁(x, sol::ShallowWaterSolution)
  op = Vector{ComplexF64}(undef, length(x))
  ∂ₓu₁!(op, x, sol)
  op
end
function ∂ₓ²u₁(x, sol::ShallowWaterSolution)
  op = Vector{ComplexF64}(undef, length(x))
  ∂ₓ²u₁!(op, x, sol)
  op
end
function ∂ₓ³u₁(x, sol::ShallowWaterSolution)
  op = Vector{ComplexF64}(undef, length(x))
  ∂ₓ³u₁!(op, x, sol)
  op
end
function u₂(x, sol::ShallowWaterSolution)
  op = Vector{ComplexF64}(undef, length(x))
  u₂!(op, x, sol)
  op
end
function ∂ₓu₂(x, sol::ShallowWaterSolution)
  op = Vector{ComplexF64}(undef, length(x))
  ∂ₓu₂!(op, x, sol)
  op
end
function ∂ₓ²u₂(x, sol::ShallowWaterSolution)
  op = Vector{ComplexF64}(undef, length(x))
  ∂ₓ²u₂!(op, x, sol)
  op
end
function ∂ₓ³u₂(x, sol::ShallowWaterSolution)
  op = Vector{ComplexF64}(undef, length(x))
  ∂ₓ³u₂!(op, x, sol)
  op
end
function solve(Ice::Ice, Fluid::Fluid, ω)
  solve(Ice, Fluid, ω, FreeClamped(), ShallowWater())
end
######################################################################################
# Function solve: Solves a frequency-domain problem
# Default option is to solve the ice-shelf problem using the FreeClamped() condition
#######################################################################################

function preallocate_matrices(::ShallowWater, ptype::Union{FreeClamped, FreeHinged})
  b = Vector{ComplexF64}(undef,1)
  c = Vector{ComplexF64}(undef,6)
  p = Vector{ComplexF64}(undef,1)
  m =  Vector{ComplexF64}(undef,6)
  a₀ = Vector{ComplexF64}(undef,1)
  ndp = preallocate_matrices(NonDimensionalProblem)
  A = Matrix{ComplexF64}(undef,7,7)
  x = Vector{ComplexF64}(undef,7)
  sw=ShallowWaterSolution(a₀, b, c, p, m,
                          ndp, A, x, ptype)
  sw,nothing
end
function solve!(cache, Ice::Ice, Fluid::Fluid, ω, ptype::Union{FreeClamped, FreeHinged}, ::ShallowWater)
  sw,nothing = cache
  non_dimensionalize!(sw.ndp, Ice, Fluid, ω)
  α = sw.ndp.α
  𝑙 = sw.ndp.𝑙
  g = sw.ndp.geo[end]
  γ = sw.ndp.γ
  𝑘 = sw.ndp.𝑘
  X = sw.ndp.X
  ω = √(α*g/𝑙)
  L = sw.ndp.geo[1]*𝑙
  H = sw.ndp.geo[2]*𝑙
  d = γ*𝑙
  Aₚ = g/(1im*ω)
  # 2) Ice-shelf
  sw.x .= 0
  sw.x[1] = α*𝑙/(H-d)
  sw.x[3] = 1-γ*α
  sw.x[end] = 1
  sw.m .= 0; PolynomialRoots.roots!(sw.m, sw.x, 1e-17, 6, true)
  sw.m = sort(sw.m, by=real, rev=true)
  m = sw.m
  # 3) Open-ocean
  k = sqrt(α*𝑙/H)
  for i=1:6
    if(i==1 || i==2)
      sw.A[i,1] = m[i]
      sw.A[i,2] = m[i]^2
      sw.A[i,3] = (ptype==FreeHinged()) ? m[i]^4 : m[i]^3
      sw.A[i,4] = m[i]^4*exp(-m[i]*L/𝑙)
      sw.A[i,5] = m[i]^5*exp(-m[i]*L/𝑙)
      sw.A[i,6] = exp(-m[i]*L/𝑙)
      sw.A[i,7] = m[i]*exp(-m[i]*L/𝑙)
    elseif(i > 2)
      sw.A[i,1] = m[i]*exp(m[i]*L/𝑙)
      sw.A[i,2] = m[i]^2*exp(m[i]*L/𝑙)
      sw.A[i,3] = (ptype==FreeHinged()) ? m[i]^4*exp(m[i]*L/𝑙) : m[i]^3*exp(m[i]*L/𝑙)
      sw.A[i,4] = m[i]^4
      sw.A[i,5] = m[i]^5
      sw.A[i,6] = 1
      sw.A[i,7] = m[i]
    end
  end
  sw.x .= 0
  sw.A[7,6] = -1
  sw.A[7,7] = 1im*k*(H/(H-d))
  sw.x[6] = Aₚ
  sw.x[7] = Aₚ*1im*k*H/(H-d)
  sw.x = transpose(sw.A)\sw.x
  for i=1:length(sw.c)
    sw.c[i] = sw.x[i]
  end
  sw.a₀[1] = sw.x[7]
  sw.BeamType = ptype
  return nothing
end
###########################################################
# u₁ is the displacement of the ice-shelf above the fluid #
###########################################################
function u₁!(op, xs, sol::ShallowWaterSolution)
  X = sol.ndp.X
  LL = sol.ndp.geo[1]
  m = sol.m
  c = sol.c
  for i in 1:length(xs)
    x = xs[i]
    if(sol.BeamType isa Union{FreeClamped, FreeHinged})
      op[i] = X*(c[1]*m[1]^2*exp.(m[1]*(x .- LL)) +
        c[2]*m[2]^2*exp.(m[2]*(x .- LL)) +
        c[3]*m[3]^2*exp.(m[3]*x) +
        c[4]*m[4]^2*exp.(m[4]*x) +
        c[5]*m[5]^2*exp.(m[5]*x) +
        c[6]*m[6]^2*exp.(m[6]*x))
    elseif(sol.BeamType isa FreeBedrock)
      X = sol.ndp.X
      xg = sol.ndp.geo[4]
      op[i] = X*(c[1]*m[1]^2*exp.(-m[1]*x) +
        c[2]*m[2]^2*exp.(-m[2]*x) +
        c[3]*m[3]^2*exp.(m[3]*x) +
        c[4]*m[4]^2*exp.(m[4]*x) +
        c[5]*m[5]^2*exp.(-m[5]*(x .- xg)) +
        c[6]*m[6]^2*exp.(-m[6]*(x .- xg)))
    end
  end
  return nothing
end

#####################################
# Freq. Domain problem with bedrock #
#####################################
function preallocate_matrices(::ShallowWater, ::FreeBedrock)
  b = Vector{ComplexF64}(undef,4)
  c = Vector{ComplexF64}(undef,6)
  p = Vector{ComplexF64}(undef,4)
  m = Vector{ComplexF64}(undef,6)
  a₀ = Vector{ComplexF64}(undef,1)
  ndp = preallocate_matrices(NonDimensionalProblem)
  A = Matrix{ComplexF64}(undef,9,9)
  x = Vector{ComplexF64}(undef,9)
  sw=ShallowWaterSolution(a₀, b, c, p, m,
                          ndp, A, x, FreeBedrock())
  cfs = Vector{ComplexF64}(undef,7)
  cfs1 = Vector{ComplexF64}(undef,5)
  sw,(cfs,cfs1)
end
function solve!(cache, Ice::Ice, Fluid::Fluid, ω, ::FreeBedrock, ::ShallowWater)
  sw,(rs,rs1) = cache
  non_dimensionalize!(sw.ndp, Ice, Fluid, ω)
  ndp = sw.ndp
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
  rs .= 0
  rs[1] = α*𝑙/(H-d)
  rs[3] = 1-γ*α
  rs[end] = 1
  #sw.m = roots(rs)
  sw.m .= 0
  PolynomialRoots.roots!(sw.m, rs, 1e-17, 6, true)
  sw.m = sort(sw.m, by=real, rev=true)
  m = sw.m
  # Solve the beam-bedrock equation
  sw.p .= 0
  rs1 .= 0
  rs1[1] = 𝑘^4 - γ*α
  rs1[end] = 1
  PolynomialRoots.roots!(sw.p, rs1, 1e-17, 4, true)
  #sw.p = roots(rs1)
  p₁ = 0; p₂ = 0
  if(real(𝑘^4 - γ*α) > 0)
    p₁ = sw.p[(real(sw.p) .> 1e-9)][1]
    p₂ = sw.p[(real(sw.p) .> 1e-9)][2]
  else
    p₁ = sw.p[abs.(real(sw.p)) .< 1e-9][1]
    p₂ = sw.p[abs.(real(sw.p)) .< 1e-9][2]
  end
  # 3) Open-ocean
  k = sqrt(α*𝑙/H)
  # Construct the linear system
  sw.A[:,1] = SA[m[1]^4, m[2]^4, m[3]^4, m[4]^4, m[5]^4*exp(m[5]*xg), m[6]^4*exp(m[6]*xg), 0, 0, 0]
  sw.A[:,2] = SA[-m[1]^5, -m[2]^5, m[3]^5, m[4]^5, -m[5]^5*exp(m[5]*xg), -m[6]^5*exp(m[6]*xg), 0, 0, 0]
  sw.A[:,3] = SA[1, 1, 1, 1, exp(m[5]*xg), exp(m[6]*xg), -1, 0, 0]
  sw.A[:,4] = SA[-m[1], -m[2], m[3], m[4], -m[5]*exp(m[5]*xg), -m[6]*exp(m[6]*xg), 1im*k*H/(H-d), 0, 0]
  sw.A[:,5] = SA[X*m[1]^2*exp(-m[1]*xg), X*m[2]^2*exp(-m[2]*xg), X*m[3]^2*exp(m[3]*xg),
                 X*m[4]^2*exp(m[4]*xg), X*m[5]^2, X*m[6]^2, 0,
                 -1, -1]
  sw.A[:,6] = SA[-X*m[1]^3*exp(-m[1]*xg), -X*m[2]^3*exp(-m[2]*xg), X*m[3]^3*exp(m[3]*xg),
                 X*m[4]^3*exp(m[4]*xg), -X*m[5]^3, -X*m[6]^3, 0,
                 p₁, p₂]
  sw.A[:,7] = SA[X*m[1]^4*exp(-m[1]*xg), X*m[2]^4*exp(-m[2]*xg), X*m[3]^4*exp(m[3]*xg),
                 X*m[4]^4*exp(m[4]*xg), X*m[5]^4, X*m[6]^4, 0,
                 -p₁^2, -p₂^2]
  sw.A[:,8] = SA[-X*m[1]^5*exp(-m[1]*xg), -X*m[2]^5*exp(-m[2]*xg), X*m[3]^5*exp(m[3]*xg),
                 X*m[4]^5*exp(m[4]*xg), -X*m[5]^5, -X*m[6]^5, 0,
                 p₁^3, p₂^3]
  sw.A[:,9] = SA[-m[1]*exp(-m[1]*xg), -m[2]*exp(-m[2]*xg), m[3]*exp(m[3]*xg), m[4]*exp(m[4]*xg),
                 -m[5], -m[6], 0, 0, 0]
  sw.x[3] = Aₚ
  sw.x[4] = 1im*Aₚ*k*H/(H-d)
  sw.x = transpose(sw.A)\sw.x
  sw.c = sw.x[1:6]
  sw.a₀ = [sw.x[7]]
  sw.b[1] = sw.x[8]; sw.b[2] = sw.x[9]
  sw.p[1] = p₁; sw.p[2] = p₂
  return nothing
end
function u₂!(op, xs, sol::ShallowWaterSolution)
  xg = sol.ndp.geo[4]
  p = sol.p
  b = sol.b
  for i in 1:length(xs)
    x = xs[i]
    op[i] = b[1]*exp.(-p[1]*(x .-xg)) + b[2]*exp.(-p[2]*(x .-xg))
  end
  return nothing
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
function ∂ₓu₁!(op, xs, sol::ShallowWaterSolution)
  X = sol.ndp.X
  LL = sol.ndp.geo[1]
  m = sol.m
  c = sol.c
  𝑙 = sol.ndp.𝑙
  for i in 1:length(xs)
    x = xs[i]
    if(sol.BeamType isa Union{FreeClamped, FreeHinged})
      op[i] = (X/𝑙)*(c[1]*m[1]^3*exp.(m[1]*(x .- LL)) +
        c[2]*m[2]^3*exp.(m[2]*(x .- LL)) +
        c[3]*m[3]^3*exp.(m[3]*x) +
        c[4]*m[4]^3*exp.(m[4]*x) +
        c[5]*m[5]^3*exp.(m[5]*x) +
        c[6]*m[6]^3*exp.(m[6]*x))
    elseif(sol.BeamType isa FreeBedrock)
      X = sol.ndp.X
      xg = sol.ndp.geo[4]
      op[i] = (X)*(-c[1]*m[1]^3*exp.(-m[1]*x) -
        c[2]*m[2]^3*exp.(-m[2]*x) +
        c[3]*m[3]^3*exp.(m[3]*x) +
        c[4]*m[4]^3*exp.(m[4]*x) -
        c[5]*m[5]^3*exp.(-m[5]*(x .- xg)) -
        c[6]*m[6]^3*exp.(-m[6]*(x .- xg)))
    end
  end
  return nothing
end
function ∂ₓu₂!(op, xs, sol::ShallowWaterSolution)
  xg = sol.ndp.geo[4]
  p = sol.p
  b = sol.b
  for i in 1:length(xs)
    x = xs[i]
    op[i] = (-p[1]*b[1]*exp.(-p[1]*(x .-xg)) - p[2]*b[2]*exp.(-p[2]*(x .-xg)))
  end
  return nothing
end

######################################
# Function to compute Bending moment #
######################################
function ∂ₓ²u₁!(op, xs, sol::ShallowWaterSolution)
  X = sol.ndp.X
  LL = sol.ndp.geo[1]
  m = sol.m
  c = sol.c
  𝑙 = sol.ndp.𝑙
  for i in 1:length(xs)
    x = xs[i]
    if(sol.BeamType isa Union{FreeClamped, FreeHinged})
      op[i] = (X)*(c[1]*m[1]^4*exp.(m[1]*(x .- LL)) +
        c[2]*m[2]^4*exp.(m[2]*(x .- LL)) +
        c[3]*m[3]^4*exp.(m[3]*x) +
        c[4]*m[4]^4*exp.(m[4]*x) +
        c[5]*m[5]^4*exp.(m[5]*x) +
        c[6]*m[6]^4*exp.(m[6]*x))
    elseif(sol.BeamType isa FreeBedrock)
      X = sol.ndp.X
      xg = sol.ndp.geo[4]
      op[i] = (X)*(c[1]*m[1]^4*exp.(-m[1]*x) +
        c[2]*m[2]^4*exp.(-m[2]*x) +
        c[3]*m[3]^4*exp.(m[3]*x) +
        c[4]*m[4]^4*exp.(m[4]*x) +
        c[5]*m[5]^4*exp.(-m[5]*(x .- xg)) +
        c[6]*m[6]^4*exp.(-m[6]*(x .- xg)))
    end
  end
  return nothing
end
function ∂ₓ²u₂!(op, xs, sol::ShallowWaterSolution)
  xg = sol.ndp.geo[4]
  p = sol.p
  b = sol.b
  for i in 1:length(xs)
    x = xs[i]
    op[i] = (p[1]^2*b[1]*exp.(-p[1]*(x .-xg)) + p[2]^2*b[2]*exp.(-p[2]*(x .-xg)))
  end
  return nothing
end


###################################
# Function to compute Shear force #
###################################
function ∂ₓ³u₁!(op, xs, sol::ShallowWaterSolution)
  X = sol.ndp.X
  LL = sol.ndp.geo[1]
  m = sol.m
  c = sol.c
  𝑙 = sol.ndp.𝑙
  for i in 1:length(xs)
    x = xs[i]
    if(sol.BeamType isa Union{FreeClamped, FreeHinged})
      op[i] = (X)*(c[1]*m[1]^5*exp.(m[1]*(x .- LL)) +
        c[2]*m[2]^5*exp.(m[2]*(x .- LL)) +
        c[3]*m[3]^5*exp.(m[3]*x) +
        c[4]*m[4]^5*exp.(m[4]*x) +
        c[5]*m[5]^5*exp.(m[5]*x) +
        c[6]*m[6]^5*exp.(m[6]*x))
    elseif(sol.BeamType isa FreeBedrock)
      X = sol.ndp.X
      xg = sol.ndp.geo[4]
      op[i] = (X)*(-c[1]*m[1]^5*exp.(-m[1]*x) -
        c[2]*m[2]^5*exp.(-m[2]*x) +
        c[3]*m[3]^5*exp.(m[3]*x) +
        c[4]*m[4]^5*exp.(m[4]*x) -
        c[5]*m[5]^5*exp.(-m[5]*(x .- xg)) -
        c[6]*m[6]^5*exp.(-m[6]*(x .- xg)))
    end
  end
  return nothing
end
function ∂ₓ³u₂!(op, xs, sol::ShallowWaterSolution)
  xg = sol.ndp.geo[4]
  p = sol.p
  b = sol.b
  for i in 1:length(xs)
    x = xs[i]
    op[i] = (-p[1]^3*b[1]*exp.(-p[1]*(x .-xg)) -p[2]^3*b[2]*exp.(-p[2]*(x .-xg)))
  end
  return nothing
end
