## Function to solve the Dispersion Equation and get one root.
function oneroot(alpha, guess)
  ans1=guess+1;
  out=guess;
  while(abs(ans1-out) > 1e-9)
    ans1=out;
    f=ans1*tanh(ans1) - alpha;
    difff=tanh(ans1)+ans1*(1/cosh(ans1))^2;
    out=ans1-f/difff;
  end
  return out;
end

## Function genalphastep to increment.
function genalphastep(start, last, step)
  length=floor((last-start)/step);
  alphastep=zeros(AbstractFloat, Int64(length+2),1);
  for m=1:Int64(length)+2
    alphastep[m]=start+(m-1)*step;
  end
  alphastep[Int64(length)+2]=abs(last);
  return alphastep;
end

## Function to perform homotopy
function homotopy(alpha, N)
  if(N==0)
    mroot=oneroot(1,1);
  else
    mroot=oneroot(1,1im*N*pi);
  end
  step=0.043;
  alphastep=genalphastep(1,abs(alpha),-step*(abs(alpha)<1.)+step*(abs(alpha)>=1.));
  for m=1:length(alphastep)
    mroot=oneroot(alphastep[m],mroot);
  end
  argstep=genalphastep(0,angle(alpha),-(π/30)*(angle(alpha)<0)+(π/30)*(angle(alpha)>=0.));
  newalphastep=zeros(ComplexF64,length(argstep),1);
  for m=2:length(argstep)
    newalphastep[m]=abs(alpha)*exp(1im*argstep[m]);
    mroot=oneroot(newalphastep[m],mroot);
  end
  return mroot;
end

## Function to solve the dispersion equation and get N roots
function dispersion_free_surface(alph, N, H)
  alpha=alph*H;
  mroots=zeros(typeof(alph),N+1,1);
  if(N==0)
    count=0;
    mroots[count+1]=homotopy(alpha,count);
  else
    count=0;
    mroots[count+1]=homotopy(alpha,count);
    count+=1;
    while(0<=1)
      mroots[count+1]=homotopy(alpha,count);
      if(abs(mroots[count+1]- (1im*count*π + alpha/(1im*count*π))) < 0.01)
        while(0<=1)
          mroots[count+1] = oneroot(alpha, 1im*count*π+alpha/(1im*count*π));
          if(abs(mroots[count+1]- (1im*count*π + alpha/(1im*count*π))) < 1e-15)
            for m=0:N-count
              mroots[count+m+1]=1im*(count+m)*π + alpha/(1im*(count+m)*π);
            end
            count=N;
            break
          end
          if(count==N)
            break
          end
          count+=1;
        end
      end
      if(count==N)
        break
      end
      count+=1;
    end
  end
  mroots=(-1im/H)*mroots;
  mroots[1]=-mroots[1];
  mroots
end
