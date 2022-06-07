function get_ref_coeff(cache, phi, NModes, k, kd, HH, dd, Γ, Ap, Amp)
  aa = zeros(Complex{Float64}, NModes+1, 1)
  dΓ = Measure(Γ,6)
  for i ∈ 1:NModes+1
    τ(x) = cos(kd[i]*(x[2]+HH))/cos(kd[i]*(HH-dd))
    aa[i] = sum(∫(τ*phi)*dΓ)
  end
  (A, M, f, g), Ref = cache
  Mt = transpose(M)
  T = inv(Mt)*A*inv(M)
  bb = T*aa + inv(Mt)*g - T*f
  c = inv(A)*(Mt*bb - g)
  copyto!(Ref, c/Ap)
end

function get_ref_modes(cache, phi, NModes, k, kd, HH, dd, Γ, Ap, Amp)
  aa = zeros(Complex{Float64}, NModes+1, 1)
  dΓ = Measure(Γ,6)
  for i ∈ 1:NModes+1
    τ(x) = cos(kd[i]*(x[2]+HH))/cos(kd[i]*(HH-dd))
    aa[i] = sum(∫(τ*phi)*dΓ)
  end
  (A, M, f, g), Ref = cache
  Mt = transpose(M)
  T = inv(Mt)*A*inv(M)
  bb = T*aa
  c = inv(A)*Mt*bb
  copyto!(Ref, c/Ap)
end
