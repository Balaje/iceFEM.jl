######################################################
# Module for interpolation of the linear system
######################################################
mutable struct FreqSpace{T1,T2}
  ω::VecOrMat{T1}
  H₀::VecOrMat{Matrix{T2}}
  f₀::VecOrMat{Matrix{T2}}
end
function FreqSpace{T}(nev::Int64, N::Int64) where T <: Float64
  ω = zeros(T, N)
  H = fill(zeros(Complex{T},nev,nev), N)
  f = fill(zeros(Complex{T},nev,1), N)
  FreqSpace{T,Complex{T}}(ω, H, f)
end
function FreqSpace{T}(nev::Int64, N::Int64) where T <: ComplexF64
  ω = zeros(T, N, N)
  H = fill(zeros(T,nev,nev), N, N)
  f = fill(zeros(T,nev,1), N, N)
  FreqSpace{T,T}(ω, H, f)
end
function InterpolateFreqDomain(C::FreqSpace, N::Tuple)
  ω0 = C.ω
  H0 = C.H₀
  f0 = C.f₀
  @assert length(N) == length(size(ω0))
  if(length(N)==1)
    return _interpolate_1d(ω0, (H0,f0), N)
  else
    return _interpolate_2d(ω0, (H0,f0), N)
  end
end
function _interpolate_1d(ω0, f, N)
  H0,f0=f
  nev=size(f0[1],1)
  itpH=interpolate((ω0,),H0,Gridded(Linear()))
  itpF=interpolate((ω0,),f0,Gridded(Linear()))
  ω1=LinRange(ω0[1],ω0[end],N[1])
  ω1, itpH.(ω1), itpF.(ω1)
end
function _interpolate_2d(ω0, f, N)
  H0,f0=f
  a=real(ω0[1]); b=real(ω0[end]);
  c=imag(ω0[1]); d=imag(ω0[end]);
  N0=size(ω0)
  nodes=(LinRange(a,b,N0[1]), LinRange(c,d,N0[2]))
  itpH=interpolate(nodes,H0,Gridded(Linear()))
  itpF=interpolate(nodes,f0,Gridded(Linear()))
  xq=LinRange(a,b,N[1])
  yq=LinRange(c,d,N[2])
  xqyq=(xq' .* ones(length(yq))) + 1im*(ones(length(xq))' .* yq)
  xqyq, itpH(xq,yq), itpF(xq,yq)
end
