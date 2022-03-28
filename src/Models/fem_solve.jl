function _get_laplace_mat_eb(f::FiniteDepthFEM, ndp::NonDimensionalProblem, beam_style, β, x₀, ω, Qϕ, χ)
  Ω = f.mesh
  Γ₃ = f.Γs[1]
  V = f.fespace
  η(x) = ω*ηₖ(x[1]-x₀, β, ndp, beam_style)[3]
  _return_matrices(Ω, Γ₃, V, Qϕ, χ, η)
end

# function _get_laplace_mat_eb(Ω, Γ, V, QΦ, χ, L, ω, ::FreeFree, offset)
#   if(μₘ==0)
#     η(x) = ω*(x[1]-0.5*L)
#     _return_matrices(Ω, Γ₃, Vh, Vh0, QΦ, χ, L, ω, η)
#   elseif(μₘ==-1)
#     η₂(x) = ω*x[1]^0
#     _return_matrices(Ω, Γ₃, Vh, Vh0, QΦ, χ, L, ω, η₂)
#   else
#     η₁(x) = ω*((cos(L*μₘ) - cosh(L*μₘ))*(sin(μₘ*(x[1]-offset)) + sinh(μₘ*(x[1]-offset)))-
#       (sin(L*μₘ) - sinh(L*μₘ))*(cos(μₘ*(x[1]-offset)) + cosh(μₘ*(x[1]-offset))))/
#       (cos(L*μₘ) - cosh(L*μₘ))
#     _return_matrices(Ω, Γ₃, Vh, Vh0, QΦ, χ, L, ω, η₁)
#   end
# end

function _return_matrices(Ω, Γ₃, V, QΦ, χ, η)
  dΩ=Measure(Ω,2); # Measure of the domains
  dΓ₃=Measure(Γ₃,6); #Interface boundary.
  a(u,v) = ∫( ∇(v)⊙∇(u) )*dΩ
  b(v) = ∫(η*v)*dΓ₃
  op=AffineFEOperator(a,b,V,V);
  K=op.op.matrix+QΦ;
  f=-(1im)*op.op.vector-χ[1,:];
  u = K\f
  ϕ = FEFunction(V, u)
  u,ϕ
end

# Function to build the reduced system
function _build_reduced_system(μ, ϕ₀, ϕₖ, ndp, Γ, V, beam_style, L₀)
  α = ndp.α
  β = 1
  γ = ndp.γ
  𝑙 = ndp.𝑙
  g = ndp.geo[end]
  ω = √(α*g/𝑙)

  dΓ=Measure(Γ,6)
  nev=length(μ)

  B=zeros(ComplexF64, nev, nev)
  K=zeros(ComplexF64, nev, nev)
  AB=zeros(ComplexF64, nev, nev)
  F=zeros(ComplexF64,nev,1)

  φ₀=FEFunction(V, ϕ₀)
  for i=1:nev
    μₘ=μ[i]
    η(x) = ηₖ(x[1]-L₀, μₘ, ndp, beam_style)[3]
    η²(x) = η(x)*η(x)
    B[i, i] = (1-γ*α)*sum(∫(η²)*dΓ)
    K[i,i]=β*μₘ^4*B[i,i]
    F[i]=(1im*ω/g)*sum(∫(η*φ₀)*dΓ)
    φₖ=FEFunction(V,ϕₖ[:,i])
    for j=1:nev
      μₘ=μ[j]
      ξ(x)=ηₖ(x[1]-L₀, μₘ, ndp, beam_style)[3]
      AB[j,i]=-(1im*ω/g)*sum(∫(ξ*φₖ)*dΓ)
    end
  end
  H=K+B+AB
  λ=H\F
  λ, K, B, AB, F
end

# function buildReducedSystem(μ, ϕ₀, ϕⱼ, α, β, γ, Γ, L, ω, V, ::FreeFree, L₀)
#   dΓ=Measure(Γ,6)
#   nev=length(μ)

#   B=zeros(ComplexF64, nev, nev)
#   K=zeros(ComplexF64, nev, nev)
#   AB=zeros(ComplexF64, nev, nev)
#   F=zeros(ComplexF64,nev,1)

#   φ₀=FEFunction(V, ϕ₀)
#   for i=1:nev
#     μₘ=μ[i]
#     if(μₘ==0)
#       η(x) = x[1]-0.5*L-L₀
#       cf = CellField(η, Γ)
#       B[i,i] = (1-γ*α)*sum(∫(cf*cf)*dΓ)
#       F[i]=(1im*ω/10)*sum(∫(η*φ₀)*dΓ)
#     elseif(μₘ==-1)
#       η₂(x) = x[1]^0
#       cf = CellField(η₂, Γ)
#       B[i,i] = (1-γ*α)*sum(∫(cf*cf)*dΓ)
#       F[i]=(1im*ω/10)*sum(∫(η₂*φ₀)*dΓ)
#     else
#       η₁(x) = ((cos(L*μₘ) - cosh(L*μₘ))*(sin(μₘ*(x[1]-L₀)) + sinh(μₘ*(x[1]-L₀))) -
#         (sin(L*μₘ) - sinh(L*μₘ))*(cos(μₘ*(x[1]-L₀)) + cosh(μₘ*(x[1]-L₀))))/
#         (cos(L*μₘ) - cosh(L*μₘ))

#       cf = CellField(η₁, Γ)
#       B[i,i]=(1-γ*α)*sum(∫(cf*cf)*dΓ)
#       F[i]=(1im*ω/10)*sum(∫(η₁*φ₀)*dΓ)
#     end

#     K[i,i]=β*μₘ^4*B[i,i];

#     φₖ=FEFunction(V,ϕⱼ[:,i])
#     for j=1:nev
#       μₘ=μ[j]
#       if(μₘ==0)
#         ξ(x) = x[1]-0.5*L-L₀
#         AB[i,j]=-(1im)*(ω/10)*sum(∫(ξ*φₖ)*dΓ)
#       elseif(μₘ==-1)
#         ξ₂(x) = x[1]^0
#         AB[i,j]=-(1im)*(ω/10)*sum(∫(ξ₂*φₖ)*dΓ)
#       else
#         ξ₁(x) = ((cos(L*μₘ) - cosh(L*μₘ))*(sin(μₘ*(x[1]-L₀)) + sinh(μₘ*(x[1]-L₀))) -
#           (sin(L*μₘ) - sinh(L*μₘ))*(cos(μₘ*(x[1]-L₀)) + cosh(μₘ*(x[1]-L₀))))/
#           (cos(L*μₘ) - cosh(L*μₘ))

#         AB[i,j]=-(1im)*(ω/10)*sum(∫(ξ₁*φₖ)*dΓ)
#       end
#     end
#   end

#   H=K+B+AB
#   λ=H\F
#   λ, K, B, AB, F
# end


# ## Get the open ocean displacement
# function ζ(X, c, ω, Ice::nonDimParameters, ::IncWave)
#   HH = Ice.HH
#   tc = Ice.tc
#   α = HH*(ω*tc)^2
#   N = length(c)
#   k = dispersionfreesurface(α, N, HH)
#   k[1] = -k[1]
#   ϕ = exp.(-k[1]*X)
#   for m=1:N
#     ϕ = ϕ + c[m]*exp.(k[m]*X)
#   end
#   ϕ
# end

# function ζ(X, c, ω, Ice::nonDimParameters, ::TransWave)
#   HH = Ice.HH
#   tc = Ice.tc
#   α = HH*(ω*tc)^2
#   N = length(c)
#   k = dispersionfreesurface(α, N, HH)
#   k[1] = -k[1]
#   ϕ = 0*X
#   for m=1:N
#     ϕ = ϕ + c[m]*exp.(-k[m]*(X.-Ice.LL))
#   end
#   ϕ
# end
