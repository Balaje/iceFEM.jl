##########################################################
# Function to Compute Resonance Frequency of the system  #
##########################################################
function computeResonanceFrequency(ice, fluid, ω₀; verbosity=0, tol=1e-8)
  count = 1
  dw = 1e-6
  fd = FiniteDepth(4)
  Δw = 1
  while (abs(Δw) > tol) && (count ≤ 50)
    s₁ = solve(ice, fluid, ω₀, FreeFree(), fd)
    s₂ = solve(ice, fluid, ω₀ + Δw, FreeFree(), fd)
    Hⁿ = s₁.K
    Hⁿ⁺¹ = s₂.K

    ΔH = (Hⁿ⁺¹ - Hⁿ)/Δw
    dwₛ = eigvals(-ΔH\Hⁿ)
    (verbosity > 0) && print(string(Δw)*"\n")
    dwₛ = sort(dwₛ,by=x->abs(x))
    Δw = dwₛ[1]
    condH = cond(Hⁿ⁺¹)
    (verbosity > 1) && print("Δw = "*string(Δw)*"\t ω₀ = "*string(ω₀)*"\t cond(H) = "*string(condH)*"...\n")
    if(!isinf(abs(Δw)))
      ω₀ = ω₀ + Δw
    end
    count+=1
  end
  (verbosity == 1) && print("Number of iterations = "*string(count)*"\n")
  ω₀
end
