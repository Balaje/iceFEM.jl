function dispersion_elastic_surface(alpha, beta, gamma, N, h)
  alpha = alpha*h
  beta = beta/h^4
  gamma = gamma/h
  del = (1 - alpha*gamma)/alpha * (beta/alpha)^(-1/5)
  H = (beta/alpha)^(-1/5)
  K = RTS_ice_roots(real(del), real(H), N+1)
  mroots = -1im*(beta/alpha)^(-1/5)*K[1:N+1]
  mroots = mroots/h
  v = mroots[3]
  mroots[3] = -mroots[1]
  mroots[1] = v
  mroots
end

function RTS_ice_roots(del, H, N)
  # Get real root
  Kgs = roots([-1, del, 0, 0, 0, 1])
  guess0 = ((real(Kgs) .> 0) .& (imag(Kgs) .== 0))'*Kgs
  K0,_ =  gen_root_ice(del, H, guess0)
  K = K0

  if(N>0)
    K = zeros(ComplexF64, N+3,1); K[1] = K0
    Del = π/H; delc0 = -Del^4; tol = 1e-8
    nc = Complex(del/delc0)^0.25
    if(del ≥ 0)
      nc = 1
    elseif nc!=round(nc)
      nc = round(nc + 0.5)
      nc = Int(abs(nc))
    end

    if nc <= N ## Get Easy imaginary roots
      w = RTS_imag_roots_ice(del, H, nc, N)
      j_notcvg = findall(w[2:length(w)]==0)
      for r=1:length(j_notcvg)
        j = j_notcvg(r)+nc; i1 = j*π; i0 = i1 - π/2
        w[j+1-nc],_ =RTS_imag_root_ice(del, H, i0, i1)
      end
      if w[1]==0
        i1 = nc*π; i0 = i1-π
        w[1],_ = RTS_imag_root_ice(del, H, i0 - tol, i1+tol)
      end
      jj = nc:N;
      K[jj.+3] = 1im*w/H
    else
      i1 = nc*π; i0 = i1-π

      K[nc.+3],_ = 1im/H*RTS_imag_root_ice(del, H, i0-tol, i1+tol)[1]
    end

    ## Get more difficult imaginary roots
    for j=1:min(N,nc-1)
      i0=(j-1)*pi; i1=i0+pi/2;
      w,_=RTS_imag_root_ice(del,H,i0,i1);
      K[j+3]=1im*w/H;
    end
    guess1 = ((real(Kgs) .> 0) .& (imag(Kgs) .> 0))⋅Kgs
    k,_ = gen_root_ice(del, H, guess1)
    tol = 1e-8
    k_cx, k_im = is_complex(k, K[nc+3], nc*π, H)
    if((k_cx .== 0) .& (k_im .== 0))
      K = complex_roots(del, H, K, nc)
    elseif k_cx .== 0
      K[3] = k_im;
      K = complex_roots(del, H, K, nc)
    else
      K[2] = k_cx;
      K[3] = -K[2]'
    end
    K = K[1:N+3]
  end
  K
end

## Given an initial guess, use Newton Raphson to find the nearest root
function gen_root_ice(del, H, guess)
  tol = 1e-8; max_reps = 50; k0 = guess
  f,df = f_df(k0, del, H)
  dk = f/df
  k = k0-dk
  reps = 0
  while (abs(dk) > tol) & (reps<max_reps)
    reps = reps+1; k0 = k
    f,df = f_df(k0, del, H)
    dk = f/df
    k = k0-dk
  end
  if reps==max_reps
    k=0
  end
  k,df
end

function f_df(K,del,H)
  Lam = K^4+ del; Lampr = 5*K^4+del; x=7.5
  if(real(K*H)<=x)
    f = Lam*K*sinh(K*H) - cosh(K*H)
    df = Lam*K*H*cosh(K*H) + (Lampr-H)*sinh(K*H)
  else
    f = Lam*K + tanh(K*H) - 1
    df = Lam*K*H + (Lampr-H)*tanh(K*H)
  end
  f,df
end
####-----END OF MAIN PROGRAM ------

function is_complex(k, vargs...)
  tol = 1e-8; y_cx = 0; y_im = 0
  k = abs(real(k)) + 1im*abs(imag(k))
  if(real(k) > tol)
    if(imag(k)>tol)
      y_cx = k
    end
  elseif length(vargs)==3
    k = 1im*sign(imag(k))*imag(k)
    Z = vargs; Kc = Z[1]; i1=Z[2]; H=Z[3]
    i0 = i1-π; w = imag(k)*H
    if((abs(k-Kc) > tol/H) & ((w-i0+tol)*(i1+tol-w) > 0))
      y_im = k
    end
  end
  y_cx, y_im
end

function complex_roots(del, H, K, nc)
  tol = 1e-8; Del=π/H; N=length(K)-3; req23 = 1
  I1 = nc*π; I0=I1-π
  k = K[nc+3]
  w = imag(k)*H
  p,dp = p_fxn(w, del, H)

  H3, del3, w3 = triple_root(nc)

  if(K[3]!=0)
    k_1 = K[3]; w_1 = imag(k_1)*H;
    p,dp_1 = p_fxn(w_1, del, H)
    req23 = 0
    if w > w_1
      W = w; dP = dp; W_1 = w_1; dP_1 = dp_1
    else
      W = w_1; dP = dp_1; W_1 = w; dP_1 = dp
    end
    if dp*dp_1 > 0
      i0 = W_1 + tol; i1 = W - tol;
      w_2 = RTS_imag_root_ice(del, H, i0, i1)[1]; K[2] = 1im*w_2/H
    elseif dP*(-1)^nc < 0
      i0 = W + tol; i1 = I1
      w_2 = RTS_imag_root_ice(del, H, i0, i1)[1]; K[2] = 1im*w_2/H
    else
      i0 = I0; i1 = W_1 - tol;
      w_2 = RTS_imag_root_ice(del, H, i0, i1)[1]; K[2] = 1im*w_2/H
    end
  elseif dp*(-1)^nc<0
    i0 = I0; i1 = w-tol; req23=0;
    _w = RTS_imag_root_ice(del, H, i0, i1)[1]; K[3] = 1im*_w/H
    i0 = w+tol; i1 = I1
    _w = RTS_imag_root_ice(del, H, i0, i1)[1]; K[2] = 1im*_w/H
  else
    if (abs(H-H3)<tol) & (abs(del-del3) < tol)
      K[2:3] = k; req23=0;
    end
  end

  w20,del20,w21,del21=double_roots(H,nc,H3,del3,w3);
  if req23==1
    if abs(del-del20)<tol
      K[2:3]=1im*w20/H; req23=0;
    elseif abs(del-del21)<tol
      K[2:3]=1im*w21/H; req23=0;
    elseif (del>del20+tol) & (del<del21-tol)
      nn=10; req23=0;
      Dw=(w21-w20)/nn; ww=w20:Dw:w21; pp=p_fxn(ww,del,H);
      if w>w21
        j=findall(pp*dp >= 0);
        while isempty(j)==1
          Dw=Dw/10; ww=w20:Dw:w21; pp=p_fxn(ww,del,H); j=findall(pp*dp >= 0);
        end
      else
        j=findall(pp*dp <= 0)
        while isempty(j)==1
          Dw=Dw/10; ww=w20:Dw:w21; pp=p_fxn(ww,del,H); j=findall(pp*dp <= 0);
        end
      end
      if length(j)==1
        if pp[j]==0
          K[2]=1im*ww[j]/H; i0=ww[j]+tol; i1=ww[j+1];
          K[3]=1im*RTS_imag_root_ice(del,H,i0,i1)/H;
        else
          i0=ww[j-1]; i1=ww[j]; K[2]=1im*RTS_imag_root_ice(del,H,i0,i1)/H;
          i0=i1; i1=ww[j+1]; K[3]=1im*RTS_imag_root_ice(del,H,i0,i1)/H;
        end
      else
        nj=length(j);
        if pp(j[1]==0 & pp[j[nj]]==0)
          K[2]=1im*ww[j[1]/H]; K[3]=1im*ww[j[nj]]/H;
        elseif pp[j[1]]==0
          K[2]=1im*ww[j[1]]/H; i0=ww[j[nj]]; i1=ww[j[nj]+1];
          K[3]=1im*RTS_imag_root_ice(del,H,i0,i1)/H;
        elseif pp[j[nj]]==0
          K[3]=1im*ww[j[nj]]/H; i0=ww[j[1]-1]; i1=ww[j[1]];
          K[2]=1im*RTS_imag_root_ice(del,H,i0,i1)/H;
        else
          i0=ww[j[1]-1]; i1=ww[j[1]];
          K[2]=1im*RTS_imag_root_ice(del,H,i0,i1)/H;
          i0=ww[j[nj]]; i1=ww[j[nj]+1];
          K[3]=1im*RTS_imag_root_ice(del,H,i0,i1)/H;
        end
      end
    end
  end

  if req23==0
    j=2:nc+3; w=sort(imag(K[j])); K[j]=1im*w;
  end

  if ((req23==1) & (nc==1) & (del<del20) & (del>-1/3/(2*H)^(2/3)))
    w=roots([-1, del*H^(2/3), 0, 1])
    k0=sqrt.(w*H^(1/3))
    k0=k0[findall(imag(k0) .> 0)];
    k,_=gen_root_ice(del,H,k0); k=is_complex(k);
    if k!=0
      K[2]=k; K[3]=-k'; req23=0;
    end
  end

  if req23==1
    _H=H .+ 0.1;
    r=roots([-1, del, 0, 0, 0, 1])
    gs=findall((real(r) .> 0) .& (imag(r) .> 0))
    gs = (length(gs) == 1) ? gs[1] : gs
    _k,df = gen_root_ice(del,_H,gs);
    _k,=is_complex(_k);
    while abs(_k)==0
      _H=_H .+ 0.1;
      r=roots([-1, del, 0, 0, 0, 1])
      gs=findall((real(r) .> 0) .& (imag(r) .> 0))
      gs = (length(gs) == 1) ? gs[1] : gs
      _k,df=gen_root_ice(del,_H,gs);
      _k,=is_complex(_k);
    end
    k=complex_seed(del,H,_k,df,_H);
    K[2]=k; K[3]=-k';
  end
  K
end

function p_fxn(w, del, H)
  p=(w.^4+del*H^4).*w.*sin(w)+H^5*cos(w);
  dp=(5*w.^4+H^4*(del-H)).*sin(w)+(w.^4+del*H^4).*w.*cos(w);
  p,dp
end

function complex_seed(del, H, vargin...)
  k = vargin[1]
  df_k = vargin[2]
  H1 = vargin[3]
  dH = -0.05
  tol = 1e-8
  if(length(vargin)==4)
    dH = vargin[4]
  end

  while H1 > H & k!=0 & abs(H1-H) ≥ tol
    H1 = H1 + dH
    Lam = k^4 + del;
    df_H = Lam*k^2*cosh(k*H) - k*sinh(k*H)
    gs = k - df_H*dH/df_k
    k,df_k = gen_root_ice(del, H, gs)
    k, = is_complex(k)
  end
  if(abs(H1-H)<tol & real(k)!=0)
    y=k
  elseif real(k) == 0
    y = complex_seed(del, H, k, df_k, H1, dH/2)
  end
  y
end

function double_roots(H, nc, vargin...)
  if(length(vargin) == 3)
    H3 = vargin[1]; del3 = vargin[2]; w3 = vargin[3]
  else
    H3, del3, w3 = triple_root(nc)
  end

  I1 = nc*π; I0 = I1-π; tol = 1e-8; I0 = maximum([I0, tol])
  jc = 2*nc - 1; Hc = ( (jc*π)^4/4 )^(0.2); wc = jc*π/2
  if(H>H3)
    w20 = Inf; w21 = Inf;
  elseif(H==H3)
    w20 = w3; w21 = w3
  else
    i0 = w3; i1 = I1; w21 = double_root(H, i0, i1)
    if H==Hc
      w20=wc
    elseif H<Hc
      i0 = I0; i1 = wc; w20 = double_root(H, i0, i1)
    else
      i0 = wc; i1 = w3; w20 = double_root(H, i0, i1)
    end
  end
  del20 = (w20 == Inf) ? NaN : -((w20/H)^4 + H*cot(w20)/w20)
  del21 = (w20 == Inf) ? NaN : -((w21/H)^4 + H*cot(w21)/w21)

  w20, del20, w21, del21
end

function double_root(H, i0, i1)
  tol = 1e-8; w = find_zero(x->q_fxn(x,H), [i0,i1], atol=tol)
  w
end

function q_fxn(w,H)
  c2 = cos(2*w); s2 = sin(2*w)
  q = 0*w - 2*H^5; j = findall(w!=0)
  q[j] = 2*w[j].^4 .*(1-c2[j]) - H^5*(1 + s2[j]/2 ./w[j])
  q
end

function triple_root(n)
  I = [(2*n-1)*π/2, n*π]
  w = find_zero(x->r_fxn(x), I, atol=1e-8)
  num = 2-2*cos(2*w) + w*sin(2*w)
  den = 2*w*cos(2*w)-sin(2*w)
  H5 = 8*w^5*num/den
  H = H5.^2; del = -( (w/H)^4 + H*cot(w)/w )
  H, del, w
end

function r_fxn(w)
  c2 = cos(2*w); s2 = sin(2*w);
  lhs = (4*w + 2*s2).*(2 - 2*c2 + w.*s2)
  rhs = (1-c2).*(2*w.*c2 - s2); y = rhs - lhs
  y
end

function RTS_imag_roots_ice(del, H, N0, N1)
  EPS=1e-9
  MAXIT=50
  x=zeros(N1+1-N0, 1);
  for j=N0:N1
    w = j*π
    i1 = w
    i0 = i1 - π/2
    if j==N0
      i0 = i1-π
    end
    its=1
    for it = 1:MAXIT
      w4 = w*w*w*w
      H4 = H*H*H*H
      p = (w4+del*H4)*w*sin(w)+H4*H*cos(w)
      dp = (5*w4 + H4*(del-H))*sin(w) + (w4 + del*H4)*w*cos(w)
      w1 = w
      w = w1 - p/dp
      if abs(w-w1) ≤ EPS
        break
      end
      its = it
    end
    if((its < MAXIT) & ((w-i0+EPS)*(i1+EPS-w) ≥ 0))
      x[j+1-N0] = w
    else
      x[j+1-N0] = 0
    end
  end
  x
end

function RTS_imag_root_ice(del, H, i0, i1)
  EPS=1e-9;
  MAXIT=100;

  H4=H*H*H*H;
  x4=i0*i0*i0*i0;
  Lam=x4+del*H4;
  fl=Lam*i0*sin(i0)+H*H4*cos(i0);

  x4=i1*i1*i1*i1;
  Lam=x4+del*H4;
  fh=Lam*i1*sin(i1)+H4*H*cos(i1);

  if (((abs(fh) > 0) & (abs(fh) > 0)) | ((abs(fl) < 0) & (abs(fh) < 0) ))
    x=0;
    exitflag=0;
    return x, exitflag
  end

  if (fl == 0)
    x = i0;
    exitflag=1;
    return x, exitflag
  end
  if (fh == 0)
    x = i1;
    exitflag=2;
    return x, exitflag
  end

  if (fl < 0.0)
    xl=i0;
    xh=i1;
  else
    xl=i1;
    xh=i0;
  end

  x=0.5*(i0+i1);
  dxold=i1-i0;
  dx=dxold;
  x4=x*x*x*x;
  Lam=x4+del*H4;
  p=Lam*x*sin(x)+H4*H*cos(x);
  dp=(5*x4+H4*(del-H))*sin(x)+Lam*x*cos(x);

  for its=1:MAXIT
    crit1=( ((x-xh)*dp-p)*((x-xl)*dp-p) >= 0 );
    crit2=( abs(2*p) > abs(dxold*dp) );
    if (crit1 | crit2)
      dxold=dx;
      dx=0.5*(xh-xl);
      x=xl+dx;
      if (x == xl)
        exitflag=3;
        return x, exitflag
      end
    else
      dxold=dx;
      dx=p/dp;
      x0=x;
      x=x0-dx;
      if x == x0
        exitflag=3;
        return x, exitflag
      end
    end
    if (abs(dx) < EPS)
      exitflag=3;
      return x, exitflag
    end
    x4=x*x*x*x;
    Lam=x4+del*H4;
    p=Lam*x*sin(x)+H4*H*cos(x);
    dp=(5*x4+H4*(del-H))*sin(x)+Lam*x*cos(x);
    if (pp < 0.0 )
      xl=x;
    else
      xh=x;
    end
  end
  x=0;
  exitflag=4;
  return x, exitflag
end
