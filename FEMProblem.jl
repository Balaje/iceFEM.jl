##############################################
# Begin file containing the FEM part
##############################################

# 1) Dispersion equation of the composite beam
function dispersion_composite_beam(beta, ndp::NonDimensionalProblem)
  α = ndp.α
  γ = ndp.γ
  𝑘 = ndp.𝑘
  xg = ndp.geo[4]

  pl = Polynomial([𝑘^4 - γ*α, 0, 0, 0, 1])
  p = roots(pl)
  p1 = 0; p2 = 0;
  if(real(𝑘^4 - γ*α) > 0)
    p1 = p[(real(p) .< 1e-9)][1]
    p2 = p[(real(p) .< 1e-9)][2]
  else
    p1 = p[abs.(real(p)) .< 1e-9][1]
    p2 = p[abs.(real(p)) .< 1e-9][2]
  end

  real( (2*beta^5 - 2*beta^5*cos(beta*xg)*cosh(beta*xg) + 2*beta*p1^2*p2^2
         - 2*beta^3*p1^2*sin(beta*xg)*sinh(beta*xg)
         - 2*beta^3*p2^2*sin(beta*xg)*sinh(beta*xg)
         + 2*beta^4*p1*cos(beta*xg)*sinh(beta*xg)
         + 2*beta^4*p1*cosh(beta*xg)*sin(beta*xg)
         + 2*beta^4*p2*cos(beta*xg)*sinh(beta*xg)
         + 2*beta^4*p2*cosh(beta*xg)*sin(beta*xg)
         - 4*beta^3*p1*p2*sin(beta*xg)*sinh(beta*xg)
         + 2*beta*p1^2*p2^2*cos(beta*xg)*cosh(beta*xg)
         - 2*beta^2*p1*p2^2*cos(beta*xg)*sinh(beta*xg)
         + 2*beta^2*p1*p2^2*cosh(beta*xg)*sin(beta*xg)
         - 2*beta^2*p1^2*p2*cos(beta*xg)*sinh(beta*xg)
         + 2*beta^2*p1^2*p2*cosh(beta*xg)*sin(beta*xg))
        /(cosh(beta*xg)*p1^2*p2^2) )
end

#######################################################
# Functions to solve the beam dispersion equations
#######################################################
function _get_roots(N, L, f::Function)
  rr=zeros(N,1)
  xbar=π/(2*L)
  count=1
  while(1>=0)
    error=1
    tol=1e-8
    r=find_zero(f,xbar)
    if(abs(r-rr[count])>tol)
      rr[count]=r
      count+=1
    end
    xbar=(count-0.5)*π/L
    if(count==N+1)
      break
    end
  end
  rr
end
function solve_eigen_eb(N, ndp::NonDimensionalProblem, ::FreeBedrock)
  xg = ndp.geo[4]
  _get_roots(N, xg, x->dispersion_composite_beam(x, ndp))
end
function solve_eigen_eb(N, ndp::NonDimensionalProblem, ::FreeClamped)
  L = ndp.geo[4]
  f(x) = cos(x*L) + 1/cosh(x*L)
  _get_roots(N, L, f)
end
function solve_eigen_eb(N, ndp::NonDimensionalProblem, ::FreeHinged)
  L = ndp.geo[4]
  f(x) = cos(x*L)*tanh(x*L) - sin(x*L)
  _get_roots(N, L, f)
end
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
  # Fix one constant and solve
  A[1,:] = [1, 0, 0, 0]
  b[1] = 1
  c = A\b
  ηg = c[1]*sinh(β*xg) + c[2]*cosh(β*xg) + c[3]*sin(β*xg) + c[4]*cos(β*xg)
  ∂ₓηg = c[1]*β*cosh(β*xg) + c[2]*β*sinh(β*xg) + c[3]*β*cos(β*xg) - c[4]*β*sin(β*xg)
  A = [1 1; p₁ p₂]
  f = [ηg, ∂ₓηg]
  b = A\f
  (c, [p₁, p₂], c[1]*sinh.(β*x) + c[2]*cosh.(β*x) + c[3]*sin.(β*x) + c[4]*cos.(β*x))
end
function ξₖ(x, β::Float64, ndp::NonDimensionalProblem, ::FreeBedrock)
  xg = ndp.geo[4]
  c, p, _ = ηₖ(xg, β, ndp)
  ηg = c[1]*sinh(β*xg) + c[2]*cosh(β*xg) + c[3]*sin(β*xg) + c[4]*cos(β*xg)
  ∂ₓηg = c[1]*β*cosh(β*xg) + c[2]*β*sinh(β*xg) + c[3]*β*cos(β*xg) - c[4]*β*sin(β*xg)
  A = [1 1; p[1] p[2]]
  f = [ηg, ∂ₓηg]
  b = A\f
  (b, b[1]*exp.(p[1]*(x .- xg)) + b[2]*exp.(p[2]*(x .- xg)))
end

function ηₖ(x, β::Float64, ndp::NonDimensionalProblem, ::FreeClamped)
  # To populate
  x
end
function ηₖ(x, β::Float64, ndp::NonDimensionalProblem, ::FreeHinged)
  # To populate
  x
end

function _plot_compare_dispersion_equation(sol₁::Union{ShallowWaterSolution, FiniteDepthSolution})
  ndp = sol₁.ndp
  xg = sol₁.ndp.geo[4]
  betaₛ = 0.001:0.001:10*π/xg
  fbetaₛ = zeros(ComplexF64,length(betaₛ),1)
  for m=1:length(betaₛ)
    fbetaₛ[m] = dispersion_composite_beam(betaₛ[m], ndp)
  end
  plt = plot(betaₛ, real.(fbetaₛ), label="Composite beam", color=:green)
  # Plot the dispersion equation of the clamped beam
  plot!(plt, betaₛ, (cos.(betaₛ*xg) + 1 ./cosh.(betaₛ*xg)), label="Clamped beam",
        line=(:dot, 1), color=:red)
  plot!(plt, betaₛ, (cos.(betaₛ*xg).*tanh.(betaₛ*xg) - sin.(betaₛ*xg)),
        label="Hinged beam", line=(:dot, 1), color=:blue)
  plot!(plt, betaₛ, 0*betaₛ, label="0", line=(:dot, 2), color=:magenta)
  title!(plt, "\$k_0 = "*string(round(𝑘,digits=4))*"\$ MPa/m")
  display(plt)
end

function _plot_invacuo_ice_eb(ndp::NonDimensionalProblem, βₛ::Vector{Float64}, modeNo::Int64)
  𝑘 = ndp.𝑘
  γ = ndp.γ
  α = ndp.α
  β = real(βₛ[modeNo])
  x₁ = 0:0.01:ndp.geo[4]
  _, _, y₁ = ηₖ(x₁, β, ndp)
  x₂ = xg:0.01:LL
  _, y₂ = ξₖ(x₂, β, ndp)
  plt = plot(x₁, real(y₁), label="\$ x < x_g\$")
  plot!(plt, x₂, real(y₂), label="\$ x > x_g\$")
  title!(plt,"\$ k^4 = "*string(round(𝑘^4,digits=4))*"\$, \$ \\gamma \\alpha = "*string(round(γ*α, digits=4))*"\$")
  plt
end
