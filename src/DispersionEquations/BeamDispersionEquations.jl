using Roots
# 1) Dispersion equation of the composite beam
function _dispersion_composite_beam(beta, ndp::NonDimensionalProblem)
  α = ndp.α
  γ = ndp.γ
  𝑘 = ndp.𝑘
  xg = ndp.geo[4]

  p = PolynomialRoots.roots([𝑘^4 - γ*α, 0, 0, 0, 1])
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
    #display(abs(r-rr[count]))
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
  _get_roots(N, xg, x->_dispersion_composite_beam(x, ndp))
end
function solve_eigen_eb(N, ndp::NonDimensionalProblem, ::FreeClamped)
  L = ndp.geo[1]
  f(x) = cos(x*L) + 1/cosh(x*L)
  _get_roots(N, L, f)
end
function solve_eigen_eb(N, ndp::NonDimensionalProblem, ::FreeHinged)
  L = ndp.geo[1]
  f(x) = cos(x*L)*tanh(x*L) - sin(x*L)
  _get_roots(N, L, f)
end
