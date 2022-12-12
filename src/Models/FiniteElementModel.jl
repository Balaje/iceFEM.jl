include("non_local_bc.jl")
include("plot_functions.jl")

########################################################
# Functions to obtain the mode shapes of the beam
########################################################
function ηₖ(x, β::Float64, ndp::NonDimensionalProblem, ::FreeBedrock)
  α = ndp.α
  γ = ndp.γ
  𝑘 = ndp.𝑘
  xg = ndp.geo[4]
  p = PolynomialRoots.roots([𝑘^4 - γ*α, 0, 0, 0, 1])
  if(real(𝑘^4 - γ*α) > 0)
    p₁ = p[(real(p) .< 1e-9)][1]
    p₂ = p[(real(p) .< 1e-9)][2]
  else
    p₁ = p[abs.(real(p)) .< 1e-9][1]
    p₂ = p[abs.(real(p)) .< 1e-9][2]
  end
  A = zeros(ComplexF64, 4, 4)
  A[1,:] = [0, 1, 0, -1];
  A[2,:] = [1, 0, -1, 0];
  A[3,:] = [β^2*sinh(β*xg) - β*cosh(β*xg)*(p₁ + p₂) + p₁*p₂*sinh(β*xg),
            β^2*cosh(β*xg) - β*sinh(β*xg)*(p₁ + p₂) + p₁*p₂*cosh(β*xg),
            p₁*p₂*sin(β*xg) - β*cos(β*xg)*(p₁ + p₂) - β^2*sin(β*xg),
            β*sin(β*xg)*(p₁ + p₂) - β^2*cos(β*xg) + p₁*p₂*cos(β*xg)]
  A[4,:] = [β^3*cosh(β*xg) - β*cosh(β*xg)*(p₁^2 + p₁*p₂ + p₂^2) + p₁*p₂*sinh(β*xg)*(p₁ + p₂),
            β^3*sinh(β*xg) - β*sinh(β*xg)*(p₁^2 + p₁*p₂ + p₂^2) + p₁*p₂*cosh(β*xg)*(p₁ + p₂),
            p₁*p₂*sin(β*xg)*(p₁ + p₂) - β*cos(β*xg)*(p₁^2 + p₁*p₂ + p₂^2) - β^3*cos(β*xg),
            β^3*sin(β*xg) + β*sin(β*xg)*(p₁^2 + p₁*p₂ + p₂^2) + p₁*p₂*cos(β*xg)*(p₁ + p₂)]
  b = zeros(ComplexF64, 4, 1)
  # Find eigenvector corresponding to λ = 0 of [A]x = λx
  ev = eigvals(A)
  ind = findfirst(abs.(ev)/maximum(abs.(ev)) .≤ 1e-8)
  c = eigvecs(A)[:,ind]

  # Find the bedrock displacement
  ηg = c[1]*sinh(β*xg) + c[2]*cosh(β*xg) + c[3]*sin(β*xg) + c[4]*cos(β*xg)
  ∂ₓηg = c[1]*β*cosh(β*xg) + c[2]*β*sinh(β*xg) + c[3]*β*cos(β*xg) - c[4]*β*sin(β*xg)
  A = [1 1; p₁ p₂]
  f = [ηg, ∂ₓηg]
  b = A\f
  (c, [p₁, p₂], (c[1]*sinh.(β*x) + c[2]*cosh.(β*x) + c[3]*sin.(β*x) + c[4]*cos.(β*x)))
end
function ξₖ(x, β::Float64, ndp::NonDimensionalProblem, ::FreeBedrock)
  xg = ndp.geo[4]
  c, p, _ = ηₖ(xg, β, ndp, FreeBedrock())
  ηg = (c[1]*sinh(β*xg) + c[2]*cosh(β*xg) + c[3]*sin(β*xg) + c[4]*cos(β*xg))
  ∂ₓηg = (c[1]*β*cosh(β*xg) + c[2]*β*sinh(β*xg) + c[3]*β*cos(β*xg) - c[4]*β*sin(β*xg))
  A = [1 1; p[1] p[2]]
  f = [ηg, ∂ₓηg]
  b = A\f
  (b, Nothing, b[1]*exp.(p[1]*(x .- xg)) + b[2]*exp.(p[2]*(x .- xg)))
end

function ηₖ(x, β::Float64, ndp::NonDimensionalProblem, ::FreeClamped)
  LL = ndp.geo[1]
  P(μ) = -(2*exp(-μ) + 2*√(2)*sin((μ+π/4)))/(exp(-2*μ) + 2*exp(-μ)*sin(μ) - 1)
  η(x,μ) = -exp.(-μ*x) - P(μ)/2*(exp.(-μ*(1 .- x)) + exp.(-μ*(1 .+ x))) - sin.(μ*x) + (1 + P(μ)*exp(-μ))*cos.(μ*x);
  (Nothing, Nothing, η(1 .- x/LL, β*LL))
end
function ηₖ(x, β::Float64, ndp::NonDimensionalProblem, ::FreeHinged)
  xg = ndp.geo[4]
  A = [0 1 0 -1;
       1 0 -1 0;
       sinh(β*xg) cosh(β*xg) sin(β*xg) cos(β*xg);
       sinh(β*xg) cosh(β*xg) -sin(β*xg) -cos(β*xg)]
  # Find eigenvector corresponding to λ = 0 of [A]x = λx
  ev = eigvals(A)
  ind = findfirst(abs.(ev)/maximum(abs.(ev)) .≤ 1e-8)
  c = eigvecs(A)[:,ind]
  # Return displacement
  (c, Nothing, (c[1]*sinh.(β*x) + c[2]*cosh.(β*x) + c[3]*sin.(β*x) + c[4]*cos.(β*x)))
end

#################################################
# Begin finite element related function.
# We use Gridap.jl to implement the method
#################################################
struct FiniteElementModel
  dim::Int64
  partition::Tuple
  nev::Int64
  NModes::Int64
end
FiniteElementModel() = FiniteElementModel(2, (20,20), 5, 3)

struct FiniteDepthFEM
  mesh::Triangulation
  Γs::Vector{BoundaryTriangulation}
  fespace::FESpace
  dim::Int64
  domain::Tuple
  partition::Tuple
  nev::Int64
end

struct FiniteElementSolution
  ϕ₀::FEFunction
  ϕₖ::Vector{FEFunction}
  λₖ::Vector{ComplexF64}
  linear_system::Tuple
  ndp::NonDimensionalProblem
  BeamType
end

function FiniteDepthFEM(ice::Ice, fluid::Fluid, ω, dim::Int64, partition::Tuple, nev, beam_type)
  @assert length(partition) == dim
  ndp = non_dimensionalize(ice, fluid, ω)
  LL = ndp.geo[1]
  HH = ndp.geo[2]
  xg = ndp.geo[4]
  γ = ndp.γ
  domain = (beam_type==FreeBedrock()) ? (0, xg, -HH, -γ) : (0, LL, -HH, -γ)
  model = CartesianDiscreteModel(domain, partition)
  model = simplexify(model)
  Ω = Triangulation(model)
  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "neumannIce", [6])
  add_tag_from_tags!(labels,"NonLocal",[7]);
  Γ₃ = BoundaryTriangulation(model, labels, tags="neumannIce");
  Γ₄ = BoundaryTriangulation(model, labels, tags="NonLocal");
  reffe = ReferenceFE(lagrangian,Float64,1);
  Vₕ=FESpace(model, reffe, conformity=:H1, vector_type=ComplexF64);
  FiniteDepthFEM(Ω, [Γ₃, Γ₄], Vₕ, dim, domain, partition, nev)
end

include("fem_solve.jl")
include("ref_coeff.jl")

function preallocate_matrices(femodel::FiniteElementModel)
  partition = femodel.partition .+ 1
  nmodes = femodel.NModes
  ndofs = partition[1]*partition[2]
  nev = femodel.nev

  m1 = spzeros(ComplexF64,ndofs,ndofs)
  stima = spzeros(ComplexF64,ndofs,ndofs)
  m2 = spzeros(ComplexF64,nmodes+1,ndofs)
  v1 = zeros(ComplexF64,ndofs)
  loadvec = zeros(ComplexF64,ndofs)

  fefunc = Vector{FEFunction}(undef, nev+1)

  H = zeros(ComplexF64,nev,nev)
  F = zeros(ComplexF64,nev)

  A = zeros(ComplexF64, nmodes+1, nmodes+1)
  M = zeros(ComplexF64, nmodes+1, nmodes+1)
  f = zeros(ComplexF64, nmodes+1)
  g = zeros(ComplexF64, nmodes+1)

  Ref = zeros(ComplexF64, nmodes+1)

  (A,M,f,g), (m1,stima,m2,v1,loadvec), (fefunc,H,F,Ref)
end

function solve!(cache, ice::Ice, fluid::Fluid, ω, ptype, femodel::FiniteElementModel; verbosity=0)
  ndp = non_dimensionalize(ice, fluid, ω)
  fem = FiniteDepthFEM(ice, fluid, ω, femodel.dim, femodel.partition, femodel.nev, ptype)
  α = ndp.α
  𝑙 = ndp.𝑙
  γ = ndp.γ
  LL = ndp.geo[1]
  HH = ndp.geo[2]
  g = fluid.g
  ω = √(α*g/𝑙)
  Aₚ = (g/(1im*ω))
  NModes = femodel.NModes

  k = dispersion_free_surface(α, NModes, HH)
  kd = dispersion_free_surface(α, NModes, HH-γ)
  μ = solve_eigen_eb(fem.nev, ndp, ptype)

  Γ₃ = fem.Γs[1]
  Γ₄ = fem.Γs[2]

  cache1, cache2, cache3 = cache

  (verbosity > 0) && print("Obtaining non-local boundary condition ...\n")
  getMQχ!((cache1,cache2), k, kd, HH, γ, NModes, Aₚ, Γ₄, fem.fespace, exp.(0*kd))
  Qϕ, K, pp, χ, f = cache2
  ϕₖʰ, H, F, Ref = cache3

  (verbosity > 0) && print("Solving Diffraction Potential ...\n")
  ϕₖʰ[1] = _get_laplace_mat_eb!((Qϕ,K,χ,f), fem, ndp, ptype, μ[1], 0, 0)
  χ=0*χ; f=0*f

  (verbosity > 0) && print("Solving Radiation Potential ... ")
  for m=1:fem.nev
    ϕₖʰ[m+1] = _get_laplace_mat_eb!((Qϕ,K,χ,f), fem, ndp, ptype, μ[m], 0, ω*𝑙)
    (verbosity > 0) && print(string(m)*"...")
  end
  (verbosity > 0) && print("\n")

  (verbosity > 0) && print("Solving the reduced system ...\n")
  λ = _build_reduced_system!(cache3, μ, ndp, Γ₃, fem.fespace, ptype, 0)
  #λ = F

  #@time _build_reduced_system!(cache3, μ, ndp, Γ₃, fem.fespace, ptype, 0)

  ϕʰ = _compute_potential(ϕₖʰ[1], ϕₖʰ[2:fem.nev+1], λ)

  get_ref_coeff((cache1,Ref), ϕʰ, NModes, k, kd, HH, γ, Γ₄, Aₚ, exp.(0*kd))

  FiniteElementSolution(ϕₖʰ[1], ϕₖʰ[2:femodel.nev+1], vec(λ),(H, F, Ref), ndp, ptype)
end

function solve(ice::Ice, fluid::Fluid, ω, ptype, femodel::FiniteElementModel; verbosity=0)
  cache = preallocate_matrices(femodel);
  solve!(cache, ice, fluid, ω, ptype, femodel; verbosity=verbosity)
end

function u₁(x, fes::FiniteElementSolution)
  λ = fes.λₖ
  ndp = fes.ndp
  nev = length(λ)
  βs = solve_eigen_eb(nev, ndp, fes.BeamType)
  U = zeros(ComplexF64, length(x), 1)
  for m=1:nev
    U = U .+ λ[m]*ηₖ(x, βs[m], ndp, fes.BeamType)[3]
  end
  U
end

function u₂(x, fes::FiniteElementSolution)
  @assert fes.BeamType == FreeBedrock()
  λ = fes.λₖ
  ndp = fes.ndp
  nev = length(λ)
  βs = solve_eigen_eb(nev, ndp, fes.BeamType)
  U = zeros(ComplexF64, length(x), 1)
  for m=1:nev
    U = U .+ λ[m]*ξₖ(x, βs[m], ndp, fes.BeamType)[3]
  end
  U
end

function _compute_potential(ϕ₀, ϕₖ, λₖ)
  ϕ = get_free_dof_values(ϕ₀)
  for i in 1:length(ϕₖ)
    ϕ = ϕ + λₖ[i]*get_free_dof_values(ϕₖ[i])
  end
  FEFunction(ϕ₀.fe_space, ϕ)
end

function ϕₕ(fes::FiniteElementSolution)
  _compute_potential(fes.ϕ₀, fes.ϕₖ, fes.λₖ)
end
