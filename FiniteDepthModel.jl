include("Input.jl")
include("DispersionFreeSurface.jl")
include("DispersionElasticSurface.jl")

struct FiniteDepth
  N::Int64
end
FiniteDepth() = FiniteDepth(3)
struct FiniteDepthSolution
  aₘ::Vector{ComplexF64}
  cₘ⁻::Vector{ComplexF64}
  cₘ⁺::Vector{ComplexF64}
  kₘ::Vector{ComplexF64}
  κₘ::Vector{ComplexF64}
  p::Vector{ComplexF64}
  b::Vector{ComplexF64}
  ndp::NonDimensionalProblem
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
  κ = dispersion_elastic_surface(α, 1., γ, N+2, HH-γ)

  function innerproduct(k, kappa, H, d)
    if(abs(k-kappa)>=1e-7)
      return ( (kappa*sin(kappa*(H-d))*cos(k*(H-d)) - k*cos(kappa*(H-d))*sin(k*(H-d)))/(kappa^2-k^2) );
    else
      return ( (H-d)/2 + sin(2*k*(H-d))/(4*k) );
    end
  end

  χ = (0:N)*(π/(HH-γ))
  λ = k; λ[1] = -k[1]

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
  if(ptype==FreeClamped)
    B4 = transpose(hcat(-κ.*exp.(-κ*LL).*tan.(κ*(HH-γ)),
                        (κ.^2).*exp.(-κ*LL).*tan.(κ*(HH-γ)),
                        -(κ.^3).*tan.(κ*(HH-γ)),
                        (κ.^4).*tan.(κ*(HH-γ))))
  else
    B4 = transpose(hcat(-κ.*exp.(-κ*LL).*tan.(κ*(HH-γ)),
                        -(κ.^3).*exp.(-κ*LL).*tan.(κ*(HH-γ)),
                        -(κ.^3).*tan.(κ*(HH-γ)),
                        (κ.^4).*tan.(κ*(HH-γ))))
  end

  B5 = transpose(hcat(-κ.*tan.(κ*(HH-γ)),
            -(κ.^2).*tan.(κ*(HH-γ)),
            -(κ.^3).*tan.(κ*(HH-γ)).*exp.(-κ*LL),
            -(κ.^4).*tan.(κ*(HH-γ)).*exp.(-κ*LL)))

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
  cₘ⁺ = cₘ[1:N+3]
  cₘ⁻ = cₘ[N+4:2N+6]
  FiniteDepthSolution(aₘ, cₘ⁺, cₘ⁻, vec(k), κ, vec(zeros(ComplexF64,2,1)),
                      vec(zeros(ComplexF64,2,1)), ndp, ptype)
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

  X = 0*x
  for m in 1:length(cₘ⁺)
    X = X + -1/(1im*ω*𝑙)*(cₘ⁻[m]*exp.(-κ[m]*x)*(-κ[m]*tan(κ[m]*(HH-γ)))
                        + cₘ⁺[m]*exp.(κ[m]*(x .-LL))*(-κ[m]*tan(κ[m]*(HH-γ))))
  end
  X
end

############################################
# Finite depth model with grounding line
############################################
