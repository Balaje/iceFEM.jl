##############################################
# Begin file containing the FEM part
##############################################
########################################################
# Functions to obtain the mode shapes of the beam
########################################################
function ηₖ(x, β::Float64, ndp::NonDimensionalProblem, ::FreeBedrock)
  α = ndp.α
  γ = ndp.γ
  𝑘 = ndp.𝑘
  xg = ndp.geo[4]
  pl = Polynomial([𝑘^4 - γ*α, 0, 0, 0, 1])
  p = roots(pl)
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
  I4 = diagm([1, 1, 1, 1])
  ev = eigvals(A, I4)
  ind = findfirst(abs.(ev)/maximum(abs.(ev)) .≤ 1e-8)
  c = eigvecs(A,I4)[:,ind]

  # Find the bedrock displacement
  ηg = c[1]*sinh(β*xg) + c[2]*cosh(β*xg) + c[3]*sin(β*xg) + c[4]*cos(β*xg)
  ∂ₓηg = c[1]*β*cosh(β*xg) + c[2]*β*sinh(β*xg) + c[3]*β*cos(β*xg) - c[4]*β*sin(β*xg)
  A = [1 1; p₁ p₂]
  f = [ηg, ∂ₓηg]
  b = A\f
  (c, [p₁, p₂], (c[1]*sinh.(β*x) + c[2]*cosh.(β*x) + c[3]*sin.(β*x) + c[4]*cos.(β*x))*(1/(c[2]+c[4])))
end
function ξₖ(x, β::Float64, ndp::NonDimensionalProblem, ::FreeBedrock)
  xg = ndp.geo[4]
  c, p, _ = ηₖ(xg, β, ndp, FreeBedrock())
  ηg = (c[1]*sinh(β*xg) + c[2]*cosh(β*xg) + c[3]*sin(β*xg) + c[4]*cos(β*xg))*(1/(c[2]+c[4]))
  ∂ₓηg = (c[1]*β*cosh(β*xg) + c[2]*β*sinh(β*xg) + c[3]*β*cos(β*xg) - c[4]*β*sin(β*xg))*(1/(c[2]+c[4]))
  A = [1 1; p[1] p[2]]
  f = [ηg, ∂ₓηg]
  b = A\f
  (b, b[1]*exp.(p[1]*(x .- xg)) + b[2]*exp.(p[2]*(x .- xg)))
end

function ηₖ(x, β::Float64, ndp::NonDimensionalProblem, ::FreeClamped)
  xg = ndp.geo[4]
  A = [0 1 0 -1;
       1 0 -1 0;
       sinh(β*xg) cosh(β*xg) sin(β*xg) cos(β*xg);
       cosh(β*xg) sinh(β*xg) cos(β*xg) -sin(β*xg)]
  # Find eigenvector corresponding to λ = 0 of [A]x = λx
  I4 = diagm([1, 1, 1, 1])
  b = zeros(ComplexF64, 4, 1)
  ev = eigvals(A, I4)
  ind = findfirst(abs.(ev)/maximum(abs.(ev)) .≤ 1e-8)
  c = eigvecs(A,I4)[:,ind]
  # Return displacement
  (c, (c[1]*sinh.(β*x) + c[2]*cosh.(β*x) + c[3]*sin.(β*x) + c[4]*cos.(β*x))*(1/(c[2]+c[4])) )
end
function ηₖ(x, β::Float64, ndp::NonDimensionalProblem, ::FreeHinged)
  xg = ndp.geo[4]
  A = [0 1 0 -1;
       1 0 -1 0;
       sinh(β*xg) cosh(β*xg) sin(β*xg) cos(β*xg);
       sinh(β*xg) cosh(β*xg) -sin(β*xg) -cos(β*xg)]
  b = zeros(ComplexF64, 4, 1)
  # Find eigenvector corresponding to λ = 0 of [A]x = λx
  I4 = diagm([1, 1, 1, 1])
  ev = eigvals(A, I4)
  ind = findfirst(abs.(ev)/maximum(abs.(ev)) .≤ 1e-8)
  c = eigvecs(A,I4)[:,ind]
  # Return displacement
  (c, (c[1]*sinh.(β*x) + c[2]*cosh.(β*x) + c[3]*sin.(β*x) + c[4]*cos.(β*x))*(1/(c[2]+c[4])))
end

function _plot_compare_dispersion_equation(sol₁::Union{ShallowWaterSolution, FiniteDepthSolution})
  ndp = sol₁.ndp
  xg = sol₁.ndp.geo[4]
  betaₛ = 0.001:0.001:10*π/xg
  fbetaₛ = zeros(ComplexF64,length(betaₛ),1)
  𝑘 = sol₁.ndp.𝑘
  for m=1:length(betaₛ)
    fbetaₛ[m] = _dispersion_composite_beam(betaₛ[m], ndp)
  end
  plt = plot(betaₛ, real.(fbetaₛ), label="Composite beam", color=:green)
  # Plot the dispersion equation of the clamped beam
  plot!(plt, betaₛ, (cos.(betaₛ*xg) + 1 ./cosh.(betaₛ*xg)), label="Clamped beam",
        line=(:dot, 1), color=:red)
  plot!(plt, betaₛ, (cos.(betaₛ*xg).*tanh.(betaₛ*xg) - sin.(betaₛ*xg)),
        label="Hinged beam", line=(:dot, 1), color=:blue)
  plot!(plt, betaₛ, 0*betaₛ, label="0", line=(:dot, 2), color=:magenta)
  title!(plt, "\$k_0 = "*string(round(𝑘,digits=4))*"\$")
  display(plt)
end

function _plot_invacuo_ice_eb(ndp::NonDimensionalProblem, modeNo::Int64)
  𝑘 = ndp.𝑘
  γ = ndp.γ
  α = ndp.α
  xg = ndp.geo[4]
  LL = ndp.geo[1]
  # Bedrock problem
  βₛ = solve_eigen_eb(10, ndp, FreeBedrock())
  β = real(βₛ[modeNo])
  x₁ = 0:0.01:ndp.geo[4]
  _, _, y₁ = ηₖ(x₁, β, ndp, FreeBedrock())
  x₂ = xg:0.01:LL

  _, y₂ = ξₖ(x₂, β, ndp, FreeBedrock())
  plt = plot(x₁, real(y₁), label="\$ x < x_g\$", color=:green, line=(:solid, 2))
  plot!(plt, x₂, real(y₂), label="\$ x > x_g\$", color=:brown, line=(:solid, 2))

  # Free Clamped beam
  βₛ = solve_eigen_eb(10, ndp, FreeClamped())
  β = real(βₛ[modeNo])
  _, y₁ = ηₖ(x₁, β, ndp, FreeClamped())
  plot!(plt, x₁, real(y₁), label="Clamped", color=:blue, line=(:dash, 1))
  plot!(plt, x₂, 0*real(x₂), label="Rigid Bed", color=:black, line=(:solid, 1))

  # Free Hinged beam
  βₛ = solve_eigen_eb(10, ndp, FreeHinged())
  β = real(βₛ[modeNo])
  _, y₁ = ηₖ(x₁, β, ndp, FreeHinged())
  plot!(plt, x₁, real(y₁), label="Hinged", color=:red, line=(:dash, 1))


  title!(plt, "Mode \$"*string(modeNo)*"\$")
  plt
end
