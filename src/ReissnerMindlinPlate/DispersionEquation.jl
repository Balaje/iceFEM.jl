## Functions to solve the Reissner Mindlin dispersion equation
# This version uses the thin-plate wave numbers as an initial guess to obtain
# roots of the Reissner Mindline version.
# (Based on the elastic dispersion code from LG Bennetts' MATLAB code (2020).
#  See DispersionIce.jl)
##

function f(z, alpha, beta, gamma, delta, zeta, H)
  #delta = 0
  A = -delta/zeta + gamma*alpha*(delta/zeta + zeta/12)
  B = (1 - gamma*alpha)*(1 - gamma*alpha*delta/12)
  C = -delta*alpha/zeta
  D = alpha*(1 - gamma*alpha*delta/12)
  # A = delta/zeta - gamma*alpha*delta/zeta
  # B = 1 - gamma*alpha
  # C = -delta*alpha/zeta
  # D = -alpha
  z.*tanh.(z*H) .- (D .+ C*z.^2)./(beta*z.^4 .+ A*z.^2 .+ B)
end

function dispersion_ice(alpha, beta, gamma, delta, zeta, N, H)
  guess_roots = -1im*dispersion_ice(alpha, beta, gamma, N, H)
  mroots = zeros(eltype(guess_roots), size(guess_roots))
  for j=1:N+1
    for m=1:length(beta)
      mroots[j,m] = nlsolve(z-> f(z, alpha, beta, gamma, delta, zeta, H),
                            [guess_roots[j,m]]; method=:newton,
                            linesearch=BackTracking(),ftol=1e-12).zero[1]
    end
  end
  mroots*1im
end
