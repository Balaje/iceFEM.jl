function solve(ice::Ice, fluid::Fluid, ω, ::FreeFree, fd::FiniteDepth, ::ReissnerMindlinIce; μ=1)
  N = fd.N
  ndp = non_dimensionalize(ice, fluid, ω, ReissnerMindlinIce(); μ=μ)
  α = ndp.α
  𝑙 = ndp.𝑙
  g = ndp.geo[5]
  γ = ndp.γ
  𝑘 = ndp.𝑘
  X = ndp.X
  ω = √(α*g/𝑙)
  LL = ndp.geo[1]
  HH = ndp.geo[2]
  d = γ*𝑙
  Aₚ = g/(1im*ω)
  # New parameters for the RM-plate
  δ = ndp.geo[end]
  ζ = γ^2

  k = dispersion_free_surface(α, N, HH)
  κ = dispersion_ice(α, 1., γ, δ, ζ, N+2, HH-γ)

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


  γαδ = (γ*α*δ/ζ + γ*α*ζ/12)
  B5 = transpose(hcat(-(κ.^3).*exp.(-κ*LL).*tan.(κ*(HH-γ)),
                      (κ.^4 .+ γαδ*κ.^2).*exp.(-κ*LL).*tan.(κ*(HH-γ)), #Free at x=L (cₘ⁻)
                      -(κ.^3).*tan.(κ*(HH-γ)),
                      (κ.^4 .+ γαδ*κ.^2).*tan.(κ*(HH-γ))) #Free at x=0 (cₘ⁻)
                 )
  B6 = transpose(hcat(-(κ.^3).*tan.(κ*(HH-γ)),
                      -(κ.^4 .+ γαδ*κ.^2).*tan.(κ*(HH-γ)), #Free at x=L (cₘ⁺)
                      -(κ.^3).*tan.(κ*(HH-γ)).*exp.(-κ*LL),
                      -(κ.^4 .+ γαδ*κ.^2).*tan.(κ*(HH-γ)).*exp.(-κ*LL)) #Free at x=0 (cₘ⁺)
                 )

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

  FiniteDepthSolution(vcat(aₘ,dₘ), cₘ⁻, cₘ⁺, vec(k), vec(κ), vec(zeros(eltype(k),2,1)),
                      vec(zeros(eltype(k),2,1)), ndp, LHS, vec(RHS), FreeFree())
end
