struct FreqSpace{T1,T2}
  ω::VecOrMat{T1}
  H::VecOrMat{Matrix{T2}}
  f::VecOrMat{Vector{T2}}
end

function FreqSpace{T}(nev::Int64, N::Int64) where T <: Float64
  ω = zeros(T, N)
  H = fill(zeros(Complex{T},nev,nev), N)
  f = fill(zeros(Complex{T},nev), N)
  FreqSpace{T,Complex{T}}(ω, H, f)
end
function FreqSpace{T}(nev::Int64, N::Int64) where T <: ComplexF64
  ω = zeros(T, N, N)
  H = fill(zeros(T,nev,nev), N, N)
  f = fill(zeros(T,nev), N, N)
  FreqSpace{T,T}(ω, H, f)
end

function interpolate_freq(ice, fluid, C::FreqSpace, N::Tuple)
end


##########################################################
# Function to Compute Resonance Frequency of the system  #
##########################################################
function computeResonanceFrequency(ice, fluid, ω₀, verbosity=0)
  count = 1
  dw = 1e-6
  tol=1e-40
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
    (verbosity > 1) && print("Δw = "*string(Δw)*"\t ω₀ = "*string(ω₀)*"\t cond(H) = "*string(condH)*"...\n")
    if(!isinf(abs(Δw)))
      ω₀ = ω₀ + Δw
    end
    count+=1
  end
  (verbosity == 1) && print("Number of iterations = "*string(count)*"\n")
  ω₀
end