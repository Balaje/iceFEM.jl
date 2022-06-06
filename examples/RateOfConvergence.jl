using iceFEM
using Plots
using Gridap

################################################
# Function for later use in the code
################################################
function ϕ(x, sol::FiniteDepthSolution)
  cₘ⁺ = sol.cₘ⁺
  cₘ⁻ = sol.cₘ⁻
  κ = sol.κₘ
  LL = sol.ndp.geo[1]
  HH = sol.ndp.geo[2]
  γ = sol.ndp.γ
  X = 0*x[1]
  for m in 1:length(cₘ⁺)
    ψₘ(z) = cos(κ[m]*(z+HH))/(cos(κ[m]*(HH-γ)))
    X = X + (cₘ⁻[m]*exp.(-κ[m]*x[1])*ψₘ(x[2])
             + cₘ⁺[m]*exp.(κ[m]*(x[1] .-LL))*ψₘ(x[2]))
  end
  X
end
function computePotential(ϕ₀, ϕₖ, λₖ)
  ϕ = get_free_dof_values(ϕ₀)
  for i in 1:length(ϕₖ)
    ϕ = ϕ + λₖ[i]*get_free_dof_values(ϕₖ[i])
  end
  FEFunction(ϕ₀.fe_space, ϕ)
end

######################################
# Begin code
######################################
ice = Ice(922.5, 2e9, 0.33, 10000, 200)
fluid = Fluid(1025, 0, 9.8, 1000, 1)

ωs = 2π*[0.04,0.02,0.01,0.005]
ns = [50,100,150,200,250,300]
hs = 1 ./ns
errs = zeros(length(ns))
plt = plot()
plot!(plt, hs, hs.^2, xaxis=:log10, yaxis=:log10, label="h²",linecolor=:black,linewidth=1,linestyle=:dash)
for (ω,clr) in zip(ωs,[:red,:green,:blue,:magenta])
  # Finite Depth Solution
  local solFD = iceFEM.solve(ice, fluid, ω, FreeClamped(), FiniteDepth(8))

  print("\n ω = "*string(ω))
  for (n,i) in zip(ns,1:length(ns))
    # Finite Element Solution
    partition = Int.((n,n/10))
    fe_model = FiniteElementModel(2, partition, 40, 8)
    solFE = iceFEM.solve(ice, fluid, ω, FreeClamped(), fe_model; verbosity=0)
    # Compute the L^2 error
    Pₕ = ϕₕ(solFE)
    ϕₑ = interpolate_everywhere(x->ϕ(x,solFD), Pₕ.fe_space)
    Ωf = get_triangulation(Pₕ)
    dΩf = Measure(Ωf, 4)
    e = abs(Pₕ-ϕₑ)
    errs[i] = √(sum(∫(e*e)*dΩf))/√(sum(∫(abs(ϕₑ)*abs(ϕₑ))*dΩf))
    print(" n = "*string(n)*"...")
  end

  plot!(plt, hs, errs, xaxis=:log10, yaxis=:log10, label="ω=2π*"*string(round(ω/2/π,digits=5)),linecolor=clr,linewidth=2)
end
