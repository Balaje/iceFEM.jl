struct FiniteDepth
  N::Int64
end
FiniteDepth() = FiniteDepth(3)
struct FiniteDepthSolution{T<:ComplexF64}
  aₘ::Vector{T}
  cₘ⁻::Vector{T}
  cₘ⁺::Vector{T}
  kₘ::Vector{T}
  κₘ::Vector{T}
  p::Vector{T}
  b::Vector{T}
  ndp::NonDimensionalProblem
  K::Matrix{T}
  f::Vector{T}
  BeamType
end

function solve(Ice::Ice, Fluid::Fluid, ω, fd::FiniteDepth)
  solve(Ice, Fluid, ω, FreeClamped(), fd)
end

function solve(Ice::Ice, Fluid::Fluid, ω, ptype::Union{FreeClamped, FreeHinged}, fd::FiniteDepth)
  N = fd.N
  ndp = non_dimensionalize(Ice, Fluid, ω)
  α = ndp.α
  𝑙 = ndp.𝑙
  g = ndp.geo[end]
  γ = ndp.γ
  𝑘 = ndp.𝑘
  X = ndp.X
  ω = √(α*g/𝑙)
  LL = ndp.geo[1]
  HH = ndp.geo[2]
  d = γ*𝑙
  Aₚ = g/(1im*ω)

  k = dispersion_free_surface(α, N, HH)
  κ = dispersion_ice(α, 1., γ, N+2, HH-γ)

  function innerproduct(k, kappa, H, d)
    if(abs(k-kappa)>=1e-7)
      return ( (kappa*sin(kappa*(H-d))*cos(k*(H-d)) - k*cos(kappa*(H-d))*sin(k*(H-d)))/(kappa^2-k^2) );
    else
      return ( (H-d)/2 + sin(2*k*(H-d))/(4*k) );
    end
  end

  χ = (0:N)*(π/(HH-γ))
  λ = k; λ[1] = k[1]

  D1 = zeros(ComplexF64, N+1, N+3)
  D2 = zeros(ComplexF64, N+1, N+1)
  D3 = zeros(ComplexF64, N+1, N+3)

  for m=1:N+1
    for n=1:N+3
      D1[m,n] = innerproduct(χ[m], κ[n], HH, γ)/(cos(χ[m]*(HH-γ))*cos(κ[n]*(HH-γ)))
      if(n ≤ N+1)
        D2[m,n] = innerproduct(χ[m], k[n], HH, γ)/(cos(χ[m]*(HH-γ))*cos(k[n]*HH))
      end
      D3[m,n] = innerproduct(k[m], κ[n], HH, γ)/(cos(k[m]*HH)*cos(κ[n]*(HH-γ)))
    end
  end

  A = diagm(vec(0.5*(cos.(k*HH).*sin.(k*HH) + k*HH)./(k.*(cos.(k*HH)).^2)))

  B1 = hcat(-D1, -D1.*transpose(repeat(exp.(-κ*LL), 1, N+1)))
  B2 = hcat(D3.*transpose(repeat(κ, 1, N+1)), D3.*transpose(repeat(-κ.*exp.(-κ*LL), 1, N+1)))
  B3 = hcat(D1.*transpose(repeat(κ.*exp.(-κ*LL), 1, N+1)), -D1.*transpose(repeat((κ), 1, N+1)))

  B4 = zeros(ComplexF64, N+1, 4)
  B5 = zeros(ComplexF64, N+1, 4)
  if(ptype==FreeClamped())
    B4 = transpose(hcat(-κ.*exp.(-κ*LL).*tan.(κ*(HH-γ)),
                        (κ.^2).*exp.(-κ*LL).*tan.(κ*(HH-γ)),
                        -(κ.^3).*tan.(κ*(HH-γ)),
                        (κ.^4).*tan.(κ*(HH-γ))))
    B5 = transpose(hcat(-κ.*tan.(κ*(HH-γ)),
                        -(κ.^2).*tan.(κ*(HH-γ)),
                        -(κ.^3).*tan.(κ*(HH-γ)).*exp.(-κ*LL),
                        -(κ.^4).*tan.(κ*(HH-γ)).*exp.(-κ*LL)))
  else
    B4 = transpose(hcat(-κ.*exp.(-κ*LL).*tan.(κ*(HH-γ)),
                        (κ.^3).*exp.(-κ*LL).*tan.(κ*(HH-γ)),
                        -(κ.^3).*tan.(κ*(HH-γ)),
                        (κ.^4).*tan.(κ*(HH-γ))))
    B5 = transpose(hcat(-κ.*tan.(κ*(HH-γ)),
                        (κ.^3).*tan.(κ*(HH-γ)),
                        -(κ.^3).*tan.(κ*(HH-γ)).*exp.(-κ*LL),
                        -(κ.^4).*tan.(κ*(HH-γ)).*exp.(-κ*LL)))
  end

  f1 = -Aₚ*D2[:,1]
  f2 = zeros(ComplexF64, N+1, 1)
  f2[1] = k[1]*Aₚ*A[1,1]

  LHS₁ = hcat(D2, B1)
  LHS₂ = hcat(diagm(vec(λ)).*A, B2)
  LHS₃ = hcat(zeros(ComplexF64,N+1,N+1), B3)
  LHS₄ = hcat(zeros(ComplexF64,4,N+1), B4, B5)

  LHS = vcat(LHS₁, LHS₂, LHS₃, LHS₄)
  RHS = vcat(f1, f2, zeros(ComplexF64, N+1, 1), zeros(ComplexF64, 4, 1))
  sol = LHS\RHS
  aₘ = sol[1:N+1]
  cₘ = sol[N+2:end]
  cₘ⁻ = cₘ[1:N+3]
  cₘ⁺ = cₘ[N+4:2N+6]
  FiniteDepthSolution(aₘ, cₘ⁻, cₘ⁺, vec(k), vec(κ), vec(zeros(ComplexF64,2,1)),
                      vec(zeros(ComplexF64,2,1)), ndp, LHS, vec(RHS), ptype)
end
############################################
# Finite depth model with grounding line
############################################
function solve(Ice::Ice, Fluid::Fluid, ω, ::FreeBedrock, fd::FiniteDepth)
  N = fd.N
  ndp = non_dimensionalize(Ice, Fluid, ω)
  α = ndp.α
  𝑙 = ndp.𝑙
  g = ndp.geo[end]
  γ = ndp.γ
  𝑘 = ndp.𝑘
  X = ndp.X
  ω = √(α*g/𝑙)
  LL = ndp.geo[1]
  HH = ndp.geo[2]
  d = γ*𝑙
  Aₚ = g/(1im*ω)
  xg = ndp.geo[4]

  k = dispersion_free_surface(α, N, HH)
  κ = dispersion_ice(α, 1., γ, N+2, HH-γ)

  function innerproduct(k, kappa, H, d)
    if(abs(k-kappa)>=1e-7)
      return ( (kappa*sin(kappa*(H-d))*cos(k*(H-d)) - k*cos(kappa*(H-d))*sin(k*(H-d)))/(kappa^2-k^2) );
    else
      return ( (H-d)/2 + sin(2*k*(H-d))/(4*k) );
    end
  end

  χ = (0:N)*(π/(HH-γ))
  λ = k; λ[1] = k[1]

  D1 = zeros(ComplexF64, N+1, N+3)
  D2 = zeros(ComplexF64, N+1, N+1)
  D3 = zeros(ComplexF64, N+1, N+3)

  for m=1:N+1
    for n=1:N+3
      D1[m,n] = innerproduct(χ[m], κ[n], HH, γ)/(cos(χ[m]*(HH-γ))*cos(κ[n]*(HH-γ)))
      if(n ≤ N+1)
        D2[m,n] = innerproduct(χ[m], k[n], HH, γ)/(cos(χ[m]*(HH-γ))*cos(k[n]*HH))
      end
      D3[m,n] = innerproduct(k[m], κ[n], HH, γ)/(cos(k[m]*HH)*cos(κ[n]*(HH-γ)))
    end
  end

  A = diagm(vec(0.5*(cos.(k*HH).*sin.(k*HH) + k*HH)./(k.*(cos.(k*HH)).^2)))

  B1 = hcat(-D1, -D1.*transpose(repeat(exp.(-κ*xg), 1, N+1)))
  B2 = hcat(D3.*transpose(repeat(κ, 1, N+1)), D3.*transpose(repeat(-κ.*exp.(-κ*xg), 1, N+1)))
  B3 = hcat(D1.*transpose(repeat(κ.*exp.(-κ*xg), 1, N+1)), -D1.*transpose(repeat((κ), 1, N+1)))


  # Solve the beam-fluid dispersion equation
  m = PolynomialRoots.roots([α/(HH-γ), 0, 1-γ*α, 0, 0, 0, 1])
  m = m[sortperm(real(m), rev=true)]
  # Solve the beam-bedrock equation
  p = PolynomialRoots.roots([𝑘^4 - γ*α, 0, 0, 0, 1])
  p₁ = 0; p₂ = 0
  if(real(𝑘^4 - γ*α) > 0)
    p₁ = p[(real(p) .> 1e-9)][1]
    p₂ = p[(real(p) .> 1e-9)][2]
  else
    p₁ = p[abs.(real(p)) .< 1e-9][1]
    p₂ = p[abs.(real(p)) .< 1e-9][2]
  end

  p₁ = -p₁
  p₂ = -p₂
  B4 = transpose(hcat((-κ.^3 - p₁*p₂*κ - (p₁+p₂)*κ.^2).*tan.(κ*(HH-γ)).*exp.(-κ*xg),
                      (κ.^4 - p₁*p₂*(p₁+p₂)*κ - (p₁^2+p₂^2+p₁*p₂)*κ.^2).*tan.(κ*(HH-γ)).*exp.(-κ*xg),
                      -(κ.^3).*tan.(κ*(HH-γ)),
                      (κ.^4).*tan.(κ*(HH-γ))))
  B5 = transpose(hcat((-κ.^3 - p₁*p₂*κ + (p₁+p₂)*κ.^2).*tan.(κ*(HH-γ)),
                      (-κ.^4 - p₁*p₂*(p₁+p₂)*κ + (p₁^2+p₂^2+p₁*p₂)*κ.^2).*tan.(κ*(HH-γ)),
                      -(κ.^3).*tan.(κ*(HH-γ)).*exp.(-κ*xg),
                      -(κ.^4).*tan.(κ*(HH-γ)).*exp.(-κ*xg)))

  f1 = -Aₚ*D2[:,1]
  f2 = zeros(ComplexF64, N+1, 1)
  f2[1] = k[1]*Aₚ*A[1,1]

  LHS₁ = hcat(D2, B1)
  LHS₂ = hcat(diagm(vec(λ)).*A, B2)
  LHS₃ = hcat(zeros(ComplexF64,N+1,N+1), B3)
  LHS₄ = hcat(zeros(ComplexF64,4,N+1), B4, B5)

  LHS = vcat(LHS₁, LHS₂, LHS₃, LHS₄)
  RHS = vcat(f1, f2, zeros(ComplexF64, N+1, 1), zeros(ComplexF64, 4, 1))
  sol = LHS\RHS
  aₘ = sol[1:N+1]
  cₘ = sol[N+2:end]
  cₘ⁻ = cₘ[1:N+3]
  cₘ⁺ = cₘ[N+4:2N+6]

  # Find the coefficients of the bedrock part
  fd = FiniteDepthSolution(aₘ, cₘ⁻, cₘ⁺, vec(k), vec(κ), vec(zeros(ComplexF64,2,1)),
                           vec(zeros(ComplexF64,2,1)), ndp, LHS, vec(RHS), FreeBedrock())
  ηg = u₁(xg, fd)
  ∂ₓηg = ∂ₓu₁(xg, fd)
  A = [1 1; p₁ p₂]
  f = [ηg, ∂ₓηg]
  b = A\f
  FiniteDepthSolution(aₘ, cₘ⁻, cₘ⁺, vec(k), vec(κ), -[p₁, p₂],
                      vec(b), ndp, LHS, vec(RHS), FreeBedrock())
end

function u₁(x, sol::FiniteDepthSolution)
  α = sol.ndp.α
  g = sol.ndp.geo[end]
  𝑙 = sol.ndp.𝑙
  ω = √(α*g/𝑙)
  cₘ⁺ = sol.cₘ⁺
  cₘ⁻ = sol.cₘ⁻
  LL = sol.ndp.geo[1]
  HH = sol.ndp.geo[2]
  γ = sol.ndp.γ
  κ = sol.κₘ
  xg = sol.ndp.geo[4]

  X = 0*x
  if(sol.BeamType isa Union{FreeClamped, FreeHinged, FreeFree})
    for m in 1:length(cₘ⁺)
      X = X + -1/(1im*ω*𝑙)*(cₘ⁻[m]*exp.(-κ[m]*x)*(-κ[m]*tan(κ[m]*(HH-γ)))
                            + cₘ⁺[m]*exp.(κ[m]*(x .-LL))*(-κ[m]*tan(κ[m]*(HH-γ))))
    end
  elseif(sol.BeamType isa FreeBedrock)
    for m in 1:length(cₘ⁺)
      X = X + -1/(1im*ω*𝑙)*(cₘ⁻[m]*exp.(-κ[m]*x)*(-κ[m]*tan(κ[m]*(HH-γ)))
                            + cₘ⁺[m]*exp.(κ[m]*(x .- xg))*(-κ[m]*tan(κ[m]*(HH-γ))))
    end
  end
  X
end
function u₂(x, sol::FiniteDepthSolution)
  xg = sol.ndp.geo[4]
  p = sol.p
  b = sol.b
  @assert length(p) == length(b)
  b[1]*exp.(-p[1]*(x .-xg)) + b[2]*exp.(-p[2]*(x .-xg))
end
######################################################
# Slope, Bending moment, Shear force
######################################################
function ∂ₓu₁(x, sol::FiniteDepthSolution)
  α = sol.ndp.α
  g = sol.ndp.geo[end]
  𝑙 = sol.ndp.𝑙
  ω = √(α*g/𝑙)
  cₘ⁺ = sol.cₘ⁺
  cₘ⁻ = sol.cₘ⁻
  LL = sol.ndp.geo[1]
  HH = sol.ndp.geo[2]
  γ = sol.ndp.γ
  κ = sol.κₘ
  xg = sol.ndp.geo[4]

  X = 0*x
  if(sol.BeamType isa Union{FreeClamped, FreeHinged})
    for m in 1:length(cₘ⁺)
      X = X + -1/(1im*ω*𝑙)*(cₘ⁻[m]*(-κ[m])*exp.(-κ[m]*x)*(-κ[m]*tan(κ[m]*(HH-γ)))
                            + cₘ⁺[m]*(κ[m])*exp.(κ[m]*(x .-LL))*(-κ[m]*tan(κ[m]*(HH-γ))))
    end
  elseif(sol.BeamType isa FreeBedrock)
    for m in 1:length(cₘ⁺)
      X = X + -1/(1im*ω*𝑙)*(cₘ⁻[m]*(-κ[m])*exp.(-κ[m]*x)*(-κ[m]*tan(κ[m]*(HH-γ)))
                            + cₘ⁺[m]*(κ[m])*exp.(κ[m]*(x .- xg))*(-κ[m]*tan(κ[m]*(HH-γ))))
    end
  end
  X
end
function ∂ₓu₂(x, sol::FiniteDepthSolution)
  xg = sol.ndp.geo[4]
  p = sol.p
  b = sol.b
  @assert length(p) == length(b)
  (-p[1]*b[1]*exp.(-p[1]*(x .-xg))
   -p[2]*b[2]*exp.(-p[2]*(x .-xg)))
end
function ∂ₓ²u₁(x, sol::FiniteDepthSolution)
  α = sol.ndp.α
  g = sol.ndp.geo[end]
  𝑙 = sol.ndp.𝑙
  ω = √(α*g/𝑙)
  cₘ⁺ = sol.cₘ⁺
  cₘ⁻ = sol.cₘ⁻
  LL = sol.ndp.geo[1]
  HH = sol.ndp.geo[2]
  γ = sol.ndp.γ
  κ = sol.κₘ
  xg = sol.ndp.geo[4]

  X = 0*x
  if(sol.BeamType isa Union{FreeClamped, FreeHinged})
    for m in 1:length(cₘ⁺)
      X = X + -1/(1im*ω*𝑙)*(cₘ⁻[m]*(-κ[m])^2*exp.(-κ[m]*x)*(-κ[m]*tan(κ[m]*(HH-γ)))
                            + cₘ⁺[m]*(κ[m])^2*exp.(κ[m]*(x .-LL))*(-κ[m]*tan(κ[m]*(HH-γ))))
    end
  elseif(sol.BeamType isa FreeBedrock)
    for m in 1:length(cₘ⁺)
      X = X + -1/(1im*ω*𝑙)*(cₘ⁻[m]*(-κ[m])^2*exp.(-κ[m]*x)*(-κ[m]*tan(κ[m]*(HH-γ)))
                            + cₘ⁺[m]*(κ[m])^2*exp.(κ[m]*(x .- xg))*(-κ[m]*tan(κ[m]*(HH-γ))))
    end
  end
  X
end
function ∂ₓ²u₂(x, sol::FiniteDepthSolution)
  xg = sol.ndp.geo[4]
  p = sol.p
  b = sol.b
  xg = sol.ndp.geo[4]
  @assert length(p) == length(b)
  ((-p[1])^2*b[1]*exp.(-p[1]*(x .-xg)) +
   (-p[2])^2*b[2]*exp.(-p[2]*(x .-xg)))
end
function ∂ₓ³u₁(x, sol::FiniteDepthSolution)
  α = sol.ndp.α
  g = sol.ndp.geo[end]
  𝑙 = sol.ndp.𝑙
  ω = √(α*g/𝑙)
  cₘ⁺ = sol.cₘ⁺
  cₘ⁻ = sol.cₘ⁻
  LL = sol.ndp.geo[1]
  HH = sol.ndp.geo[2]
  γ = sol.ndp.γ
  κ = sol.κₘ
  xg = sol.ndp.geo[4]

  X = 0*x
  if(sol.BeamType isa Union{FreeClamped, FreeHinged})
    for m in 1:length(cₘ⁺)
      X = X + -1/(1im*ω*𝑙)*(cₘ⁻[m]*(-κ[m])^3*exp.(-κ[m]*x)*(-κ[m]*tan(κ[m]*(HH-γ)))
                            + cₘ⁺[m]*(κ[m])^3*exp.(κ[m]*(x .-LL))*(-κ[m]*tan(κ[m]*(HH-γ))))
    end
  elseif(sol.BeamType isa FreeBedrock)
    for m in 1:length(cₘ⁺)
      X = X + -1/(1im*ω*𝑙)*(cₘ⁻[m]*(-κ[m])^3*exp.(-κ[m]*x)*(-κ[m]*tan(κ[m]*(HH-γ)))
                            + cₘ⁺[m]*(κ[m])^3*exp.(κ[m]*(x .- xg))*(-κ[m]*tan(κ[m]*(HH-γ))))
    end
  end
  X
end
function ∂ₓ³u₂(x, sol::FiniteDepthSolution)
  xg = sol.ndp.geo[4]
  p = sol.p
  b = sol.b
  xg = sol.ndp.geo[4]
  @assert length(p) == length(b)
  ((-p[1])^3*b[1]*exp.(-p[1]*(x .-xg)) +
   (-p[2])^3*b[2]*exp.(-p[2]*(x .-xg)))
end

### Add a new interface to implement the iceberg vibration problem
## NOTE: This file needs to be split into modules.
function solve(ice::Ice, fluid::Fluid, ω, ::FreeFree, fd::FiniteDepth)
  N = fd.N
  ndp = non_dimensionalize(ice, fluid, ω)
  α = ndp.α
  𝑙 = ndp.𝑙
  g = ndp.geo[end]
  γ = ndp.γ
  𝑘 = ndp.𝑘
  X = ndp.X
  ω = √(α*g/𝑙)
  LL = ndp.geo[1]
  HH = ndp.geo[2]
  d = γ*𝑙
  Aₚ = g/(1im*ω)

  k = dispersion_free_surface(α, N, HH)
  κ = dispersion_ice(α, 1., γ, N+2, HH-γ)

  function innerproduct(k, kappa, H, d)
    if(abs(k-kappa)>=1e-7)
      return ( (kappa*sin(kappa*(H-d))*cos(k*(H-d)) - k*cos(kappa*(H-d))*sin(k*(H-d)))/(kappa^2-k^2) );
    else
      return ( (H-d)/2 + sin(2*k*(H-d))/(4*k) );
    end
  end

  χ = (0:N)*(π/(HH-γ))
  λ = k; λ[1] = k[1]

  D1 = zeros(ComplexF64, N+1, N+3)
  D2 = zeros(ComplexF64, N+1, N+1)
  D3 = zeros(ComplexF64, N+1, N+3)

  for m=1:N+1
    for n=1:N+3
      D1[m,n] = innerproduct(χ[m], κ[n], HH, γ)/(cos(χ[m]*(HH-γ))*cos(κ[n]*(HH-γ)))
      if(n ≤ N+1)
        D2[m,n] = innerproduct(χ[m], k[n], HH, γ)/(cos(χ[m]*(HH-γ))*cos(k[n]*HH))
      end
      D3[m,n] = innerproduct(k[m], κ[n], HH, γ)/(cos(k[m]*HH)*cos(κ[n]*(HH-γ)))
    end
  end

  A = diagm(vec(0.5*(cos.(k*HH).*sin.(k*HH) + k*HH)./(k.*(cos.(k*HH)).^2)))

  B1 = hcat(-D1, -D1.*transpose(repeat(exp.(-κ*LL), 1, N+1))) # Match potentials at x=0
  B2 = hcat(D3.*transpose(repeat(κ, 1, N+1)), D3.*transpose(repeat(-κ.*exp.(-κ*LL), 1, N+1))) # Match vel. at x=0
  B3 = hcat(-D1.*transpose(repeat(exp.(-κ*LL), 1, N+1)), -D1) # Match potentials at x=L
  B4 = hcat(D3.*transpose(repeat(-κ.*exp.(-κ*LL), 1, N+1)), D3.*transpose(repeat(κ, 1, N+1))) # Match vel. at x=L


  B5 = transpose(hcat(-(κ.^3).*exp.(-κ*LL).*tan.(κ*(HH-γ)),
                      (κ.^4).*exp.(-κ*LL).*tan.(κ*(HH-γ)), #Free at x=L (cₘ⁻)
                      -(κ.^3).*tan.(κ*(HH-γ)),
                      (κ.^4).*tan.(κ*(HH-γ)))) #Free at x=0 (cₘ⁻)
  B6 = transpose(hcat(-(κ.^3).*tan.(κ*(HH-γ)),
                      -(κ.^4).*tan.(κ*(HH-γ)), #Free at x=L (cₘ⁺)
                      -(κ.^3).*tan.(κ*(HH-γ)).*exp.(-κ*LL),
                      -(κ.^4).*tan.(κ*(HH-γ)).*exp.(-κ*LL))) #Free at x=L (cₘ⁺)

  Z = zeros(ComplexF64, N+1, N+1)

  f1 = -Aₚ*D2[:,1]
  f2 = zeros(ComplexF64, N+1, 1)
  f2[1] = k[1]*Aₚ*A[1,1]

  LHS₁ = hcat(D2, B1, Z)
  LHS₂ = hcat(diagm(vec(λ)).*A, B2, Z)
  LHS₃ = hcat(Z, B3, D2)
  LHS₄ = hcat(Z, B4, diagm(vec(λ)).*A)
  LHS₅ = hcat(zeros(ComplexF64,4,N+1), B5, B6, zeros(ComplexF64,4,N+1))

  LHS = vcat(LHS₁, LHS₂, LHS₃, LHS₄, LHS₅)
  RHS = vcat(f1, f2, zeros(ComplexF64, N+1, 1),
             zeros(ComplexF64, N+1, 1), zeros(ComplexF64, 4, 1))
  sol = LHS\RHS
  aₘ = sol[1:N+1]
  cₘ = sol[N+2:3N+7]
  cₘ⁻ = cₘ[1:N+3]
  cₘ⁺ = cₘ[N+4:2N+6]
  dₘ = sol[3N+8:4N+8]

  FiniteDepthSolution(vcat(aₘ,dₘ), cₘ⁻, cₘ⁺, vec(k), vec(κ), vec(zeros(ComplexF64,2,1)),
                      vec(zeros(ComplexF64,2,1)), ndp, LHS, vec(RHS), FreeFree())
end
