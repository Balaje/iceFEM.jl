### Function to solve the ice-dispersion equation.
# Converted to Julia from LG Bennetts' MATLAB code (2020)
function oneroot(alpha, beta, gamma, guess)
  ans1 = guess+1
  out = guess
  while abs(ans1-out) > 1e-12
    ans1 = out
    out = ans1 - f(ans1, alpha, beta, gamma)/difff(ans1, alpha, beta, gamma)
  end
  out
end

function f(z, alpha, beta, gamma)
  z*tanh(z) - alpha/(beta*z^4 + (1-gamma*alpha))
end
function difff(z, alpha, beta, gamma)
  tanh(z) + z*sech(z)^2 + alpha*beta*4*z^3/(beta*z^4 + (1-alpha*gamma))^2;
end



## Homotopy: Calculates the Nth root using the homotopy method

function homotopy(alpha, beta, gamma, N)
  if(N==0)
    mroot = oneroot(1, beta, gamma, 1)
  elseif(N==1)
    arr = collect(0:0.01*(2*(beta > 1)-1):log(beta))
    betastep = exp.([arr; log(beta)])
    mroot = oneroot(1, 1, gamma, 0.6066 + 0.8761im)
    for k=2:length(betastep)
      mroot = oneroot(1, betastep[k], gamma, mroot)
    end
  elseif(N==2)
    arr = collect(0:0.01*(2*(beta > 1)-1):log(beta))
    betastep = exp.([arr; log(beta)])
    mroot = oneroot(1, 1, gamma, -0.6066 + 0.8761im)
    for k=2:length(betastep)
      mroot = oneroot(1, betastep[k], gamma, mroot)
    end
  else
    mroot = oneroot(1, beta, gamma, 1im*(N-2)*π)
  end

  arr = collect(0:0.01*(2*(abs(alpha) > 1)-1):log(abs(alpha)))
  arr1 = [arr; log(abs(alpha))]
  alphastep = exp.(arr1)

  for k=2:length(alphastep)
    mroot = oneroot(alphastep[k],beta,gamma,mroot)
  end

  if angle(alpha) > 0
    arr = collect(0:π/30:angle(alpha))
    alphastep = abs(alpha)*exp.(1im*[arr; angle(alpha)])
  else
    arr = collect(0:-π/30:angle(alpha))
    alphastep = abs(alpha)*exp.(1im*[arr; angle(alpha)])
  end

  for k=2:length(alphastep)
    mroot = oneroot(alphastep[k],beta,gamma,mroot);
  end

  mroot
end

function dispersion_ice(alpha, beta, gamma, N, H)
  alpha = alpha*H
  beta = beta/H^4
  gamma = gamma/H
  mroots = zeros(ComplexF64,N+1, length(beta))
  for j=1:N+1
    mroots[j,1] = homotopy(alpha[1], beta[1], gamma[1], j-1)
    for p=2:length(beta)
      mroots[j,p] = oneroot(alpha[p], beta[p], gamma[p], mroots[j,p-1])
    end
  end
  -(1im)*mroots/H
end
