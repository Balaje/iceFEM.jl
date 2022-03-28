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
