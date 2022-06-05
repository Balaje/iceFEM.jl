# Function to get the associated matrices on the semi--infinite domain
function innerproduct(k, kappa, H, d)
  if(abs(k-kappa)>=1e-7)
    return ( (kappa*sin(kappa*(H-d))*cos(k*(H-d)) - k*cos(kappa*(H-d))*sin(k*(H-d)))/(kappa^2-k^2) );
  else
    return ( (H-d)/2 + sin(2*k*(H-d))/(4*k) );
  end
end

# Get Matrices in the open ocean.
function getMAT(k, kd, H, d, NModes, Ap, Amp)
  A = zeros(Complex{Float64},NModes+1,NModes+1);
  M = zeros(Complex{Float64},NModes+1,NModes+1);
  f = zeros(Complex{Float64},NModes+1,1);
  g = zeros(Complex{Float64},NModes+1,1);
  for i=1:NModes+1
    A[i,i]=0.5*Amp[i]*(cos(k[i]*H)*sin(k[i]*H) + k[i]*H)/(cos(k[i]*H))^2
    f[i]=Ap*innerproduct(k[1], kd[i], H, d)/(cos(kd[i]*(H-d))*cos(k[1]*H));
    for j=1:NModes+1
      M[i,j]=innerproduct(k[j], kd[i], H, d)/(cos(kd[i]*(H-d))*cos(k[j]*H));
    end
  end
  g[1]=-Ap*A[1,1]
  A, M, f, g
end

# Get the non-local matrix on the boundary
function getMQχ!(Qϕ, χ, pp, k, kd, H, d, NModes, Ap, Γ, V, Amp)
  A,M,f,g=getMAT(k,kd,H,d,NModes,Ap,Amp);
  dΓ=Measure(Γ,8);
  for m=1:NModes+1
    τ(x)=cos(kd[m]*(x[2]+H))/cos(kd[m]*(H-d));
    a(u,v)=∫(u*v)*dΓ; #Dummy bilinear form ?
    b(v)=∫(τ*v)*dΓ;
    op=AffineFEOperator(a,b,V,V);
    pp[m,:]=op.op.vector;
  end
  # Get the matrix corresponding to Qϕ
  Mt=transpose(M);
  T=sparse(inv(Mt)*A*inv(M));
  copyto!(Qϕ,transpose(pp)*T*pp);
  c=inv(Mt)*g-T*f;
  copyto!(χ,transpose(c)*pp)
end
