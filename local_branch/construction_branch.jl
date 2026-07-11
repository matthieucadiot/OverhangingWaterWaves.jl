
### code used to produce the approximate branch
ϵ = [0.5 * mid(ϵ0) + 0.5*cospi(2*k / N_fft_)*mid(ϵ0) for k ∈ 0:npts_-1] 

T = Sequence{CartesianProduct{Tuple{ParameterSpace, CosFourier{Int64}}}, Vector{Float64}}
Uc_fft = Vector{T}(undef, npts_)

U0  = Sequence(ParameterSpace()×CosFourier(K0,1),big.(ones(K0+2)))
U0 = Newton_method_bar(U0, mid(m̄) , mid(γ), mid(g), mid(h) , 50, 1e-15)
Uc_fft[npts_] = U0

for n=npts_-1:-1:1
     U = big.(Uc_fft[n+1])
     Uc_fft[n] = Newton_method(U, ϵ[n], mid(m̄) , mid(γ), mid(g), mid(h) , 50, 1e-16)
end

U_fft= [Uc_fft ; reverse(Uc_fft)[2:npts_-1]] 
U_cheb = Vector{Sequence{Chebyshev,Vector{Interval{BigFloat}}}}(undef, K0+2)
for n=1:K0+2
    U_cheb[n] = interval.(mid.(rifft!(complex.([coefficients(U_fft[j])[n] for j = 1:2*(npts_-1)]), Chebyshev(N))))  #computation of the Chebyshev coefficients
end
U_fft_big = cheb2grid(U_cheb, N_fft_) ### this ensures that U_fft is the grid representation of a Chebyshev polynomial of order N, represented on a grid with N_fft_ points. The Cheb coefficients are given in U_cheb 
U_fft_big = [Sequence(ParameterSpace()×CosFourier(K0,interval(1)), U_fft_big[i]) for i = 1:2*(npts_-1)]   




#### code used to save the data 
jldsave("U_fft.jld2" ; U_fft_big)

  m = -I"0.1"
 Q̄ = interval(2)*g*h + (m/h - γ*h*I"0.5")^2
#### flat solution at the bifurcation point
 ā = Sequence(Sa, interval.(zeros(K0+1))) ; ā[0] = h
 Ū = Sequence(S,[Q̄ ; ā[0:K0]])

 ## we save the end point of the branch. We will use it to continue the branch further
 U = Ū + ϵ0^2*op_P(ϵ0,S)*U_fft_big[1]
 U = interval.(Float64.(inf.(U),RoundDown),Float64.(sup.(U),RoundUp) )

jldsave("m_end.jld2" ; m)
jldsave("U_end.jld2" ; U)