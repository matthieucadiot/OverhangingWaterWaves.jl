function cos2exp(a)
    N = size(a)[1]-1
    return Sequence(Fourier(N,interval(1)), vec([reverse(a[1:N]); a[0:N]]))
end

function sin2exp(a)
    N = size(a)[1]
    return Sequence(Fourier(N,interval(1)), vec([im*reverse(a[1:N]); interval(0) ; -im*a[1:N]]))
end


function Dy_op_inv(Sc,d::Interval{Float64})
    K0 = order(Sc)
    dn = interval.(1:K0)
    dn = [d ; interval(1) ./ dn]
    return  LinearOperator(Sc, Sc, Diagonal( dn ))
end 





function Dy_op_1D(Sc,d::Interval{Float64})
    K = order(Sc)
    dn = interval.(1:K)
    dn = [interval(1)/d ; dn.*coth.(d*dn)] 

    return LinearOperator(Sc, Sc, Diagonal( dn ))
end 


function Dy_op_1D(Sc,d::Float64)
    K = order(Sc)
    dn = (1:K)
    dn = [(1)/d ; dn.*coth.(d*dn)] 

    return LinearOperator(Sc, Sc, Diagonal( dn ))
end 



function D0_op(Sc::CartesianProduct{Tuple{ParameterSpace, CosFourier{Interval{Float64}}}})
    K = order(Sc[2])
    dn = [interval(1) ; interval(1) ;interval(2)*interval.(2:K+1)]
    return LinearOperator(Sc, Sc, Diagonal( dn ))
end 

function D0_op(Sc::CosFourier{Interval{Float64}})
    K = order(Sc)
    dn = [interval(1) ;interval(2)*interval.(2:K+1)]
    return LinearOperator(Sc, Sc, Diagonal( dn ))
end 

function D0_inv_op(Sc::CartesianProduct{Tuple{ParameterSpace, CosFourier{Interval{Float64}}}})
    K = order(Sc[2])
    dn = [interval(1) ; interval(1) ;interval(0.5) ./ interval.(2:K+1)]
    return LinearOperator(Sc, Sc, Diagonal( dn ))
end 

function D0_inv_op(Sc::CosFourier{Interval{Float64}})
    K = order(Sc)
    dn = [interval(1) ;interval(0.5) ./ interval.(2:K+1)]
    return LinearOperator(Sc, Sc, Diagonal( dn ))
end 




function D0_fct(U::Sequence{CartesianProduct{Tuple{ParameterSpace, CosFourier{Interval{Float64}}}}, Vector{Interval{Float64}}})
    Sc = space(U) ; K = order(Sc[2])
    dn = [interval(1) ; interval(1) ; interval(2)*interval.(2:K+1)]
    return dn.*U
end 


function D0_fct(U::Sequence{CosFourier{Interval{Float64}}, Vector{Interval{Float64}}})
    Sc = space(U) ; K = order(Sc)
    dn = [interval(1) ; interval(2)*interval.(2:K+1)]
    return dn.*U
end 


function  comp_lamb(N)
    
    f = interval(0);
    for n = 1:N
        f = f + interval(1)/interval(n)*tan(interval(π)*interval(n)/(interval(2*N+1)));
    end
    f = interval(1)/interval(2*N+1) + interval(2)/interval(π)*f;
end




function Dy_fct_1D(U,d::Interval{Float64})
    Sc = space(U)
    K = order(Sc)
    dn = interval.(1:K)
    dn = [interval(1)/d ; dn.*coth.(d*dn)]
    return dn.*U
end


#### extended precision version, used only for the residual in the Y bound
function Dy_fct_1D(U,d::Interval{BigFloat})
    Sc = space(U)
    K = order(Sc)
    dn = interval.(big.(1:K))
    dn = [interval(1)/d ; dn.*coth.(d*dn)]
    return dn.*U
end



function Dy_fct_1D(U,d::AbstractFloat)
    Sc = space(U)
    K = order(Sc)
    dn = (1:K)
    dn = [(1)/d ; dn.*coth.(d*dn)]
    return dn.*U
end 




function F_grid(m,h,γ,g,U,S,npts_)
    Tel = eltype(coefficients(U[1]))   #### Interval{Float64} or, for the extended-precision residual, Interval{BigFloat}
    T = CartesianProduct{Tuple{ParameterSpace, CosFourier{Interval{Float64}}}}
    F_grid = Vector{Sequence{T,Vector{Tel}}}(undef, npts_)
    a = [component(U[n],2) for n = 1:npts_]
    d = h
    Q = [component(U[n],1)[1] for n = 1:npts_]
    e0 = Sequence(space(a[1]), zeros(Tel, size(a[1])[1]))
    e0[0] = one(Tel)
    for n = 1:npts_
        F_grid0 = Sequence(S, zeros(Tel, order(S[2])+2))
        component(F_grid0,1)[1] .= a[n][0] - h
        component(F_grid0,2) .= project( (m[n]*Dy_fct_1D(e0,d) + interval(0.5)*γ*Dy_fct_1D(a[n]^2,d) - γ*a[n]*(Dy_fct_1D(a[n],d)))^2 - (Q[n] - interval(2)*g*a[n])*((Derivative(1)*a[n])^2 + Dy_fct_1D(a[n],d)^2), S[2])
        F_grid[n] = F_grid0
    end
    return F_grid
end



function DF_grid(m,h,γ,g,U,S,npts_)
    T = CartesianProduct{Tuple{ParameterSpace, CosFourier{Interval{Float64}}}}
    DF_grid = Vector{LinearOperator{T,T,Matrix{Interval{Float64}}}}(undef, npts_)
    a = [component(U[n],2) for n = 1:npts_]
     d = h
    Q = [component(U[n],1)[1] for n = 1:npts_]
    for n = 1:npts_
        Ma = project(Multiplication(a[n]), S, S)
        Dy = Dy_op_1D(S,d)
        C = γ*(Dy*Ma - Ma*Dy - project(Multiplication(Dy_fct_1D(a[n],d)), S, S))

        K = order(S)
        Dx = project(Derivative(1), S, SinFourier(K,interval(1)))

        e0 = Sequence(space(a[n]), interval.(zeros(size(a[n])[1])))
        e0[0] = interval(1)

        v1 = interval(2)*(m[n]*Dy*e0 + interval(0.5)γ*Dy_fct_1D(a[n]^2,d) - γ*a[n]*Dy_fct_1D(a[n],d))
        v2 = interval(2)*g*((Derivative(1)*a[n])^2 + (Dy_fct_1D(a[n],d))^2)
        v3 = -interval(2)*(Q[n] - interval(2)*g*(a[n]))*(Derivative(1)*a[n])
        v4 = -interval(2)*(Q[n] - interval(2)*g*(a[n]))*(Dy_fct_1D(a[n],d)) 

        DF0 = LinearOperator(ParameterSpace()×S, ParameterSpace()×S, interval.(zeros(K+2, K+2)))
        e0 = interval.(zeros(1,K+1)) ; e0[1] = interval(1.0)
        component(DF0,1,2).= e0
        component(DF0,2,1) .= - project((Dx*a[n])^2 + (Dy_fct_1D(a[n],d))^2, S)
        component(DF0,2,2) .= project(Multiplication(v1), S, S)*C + project(Multiplication(v2), S, S) + project(Multiplication(v3), SinFourier(K,1), S)*Dx + project(Multiplication(v4), S, S)*Dy
        DF_grid[n] = DF0
    end
    return DF_grid
end




function F(m,γ,k,g,h,U)

    Q = component(U,1)[1]
    a = component(U,2)
    d = k*h
    T = eltype(coefficients(U))

    e0 = Sequence(space(a), zeros(T, size(a)[1]))
    e0[0] = one(T)
    F = Sequence(space(U), zeros(T, size(U)[1]))

    component(F,1)[1] .= a[0] - h
    component(F,2) .= project( (m*Dy_fct_1D(e0,d) + γ/2*Dy_fct_1D(a^2,d) - γ*a*(Dy_fct_1D(a,d)))^2 - (Q-2*g*a)*((Derivative(1)*a)^2 + (Dy_fct_1D(a,d))^2) , space(a))
    return F
end



function DF(m,γ,k,g,h,U)

    d = k*h

      Q = component(U,1)[1]
    a = component(U,2) ; Sa = space(a)
    T = eltype(coefficients(U))
    N = size(a)[1]-1
    n = 1:N

Dy = LinearOperator(Sa, Sa, Diagonal( vec([1/d ; n.*coth.(d*n)]) ))
Dx = project(Derivative(1), Sa, SinFourier(N,1))

    N = order(space(a))
    Ma = project(Multiplication(a), space(a), space(a))
    C = γ*(Dy*Ma - Ma*Dy - project(Multiplication(Dy_fct_1D(a,d)), space(a), space(a)))

    e0 = Sequence(space(a), zeros(T, size(a)[1]))
    e0[0] = one(T)

    v1 = 2*(m*Dy*e0 + γ/2*Dy_fct_1D(a^2,d) - γ*a*(Dy_fct_1D(a,d)))
    v2 = 2*g*((Dx*a)^2 + (Dy_fct_1D(a,d))^2)
    v3 = -2*(Q - 2*g*(a))*(Dx*a)
    v4 = -2*(Q - 2*g*(a))*(Dy_fct_1D(a,d)) 

    DF = LinearOperator(space(U), space(U), zeros(T, size(U)[1], size(U)[1]))
    e0 = zeros(T,1,N+1) ; e0[1] = one(T)
    component(DF,1,2).= e0
    component(DF,2,2).= project(Multiplication(v1), space(a), space(a))*C + project(Multiplication(v2), space(a), space(a)) + project(Multiplication(v3), SinFourier(N,1), space(a))*Dx + project(Multiplication(v4), space(a), space(a))*Dy
    component(DF,2,1) .= - project((Dx*a)^2 + (Dy_fct_1D(a,d))^2, space(a)) 
    return DF
end




function Newton_method(m,γ,k,g,h,U, max_iter, tol)

    F_U = F(m,γ,k,g,h,U)
    DF_U = DF(m,γ,k,g,h,U)
    if norm(F_U) < tol
            println("Converged in 0 iterations.")
            return U
    end
    for iter = 1:max_iter
        
        U = U - DF_U \ F_U

        F_U = F(m,γ,k,g,h,U)
        DF_U = DF(m,γ,k,g,h,U)

        nf = norm(F_U,1)
        # println("Iteration $iter: ||F(U)|| = $nf")

        if nf < tol
            # println("Converged in $iter iterations.")
            return U
        end
    end
    # println("Did not converge within $max_iter iterations.")
    return U
end
 


function norm_seq_grid(U::Vector{Sequence{CosFourier{Interval{Float64}}, Vector{Interval{Float64}}}},N_)
    #### computation of an upper bound for the supremum of a seuquence norm. The upper bound is achieved using the minimum of 2 upper bounds.

    ### For the first one, we compute the ℓ1 norm of the chebyshev coefficients of U, where U is the representation on the grid. We first compute the chebyshev coefficients rigorously up to order N_ and then compute the norm.

    #### for the second, we use the Lebesgue constant, and take the supremum on the grid.
    S = space(U[1])
    K = order(S)
    npts_= length(U)
    U_cheb = Vector{Sequence{Chebyshev,Vector{Interval{Float64}}}}(undef, K+1)
    for k = 0:K 
         U_cheb[k+1] = rifft!(complex.([U[j][k] for j = 1:npts_]), Chebyshev(N_))  #computation of the Chebyshev coefficients
    end
    weight = interval.(2*ones(length(U_cheb))) ; weight[1] = interval(1)

    ΛN = comp_lamb(N_)
    return minimum([norm(norm.(weight.*U_cheb,1),1) ΛN*norm(norm.(U,1),Inf)])
end





function norm_seq_grid(U::Vector{Sequence{Fourier{Interval{Float64}}, Vector{Interval{Float64}}}},N_)
    #### computation of an upper bound for the supremum of a seuquence norm. The upper bound is achieved using the minimum of 2 upper bounds.

    ### For the first one, we compute the ℓ1 norm of the chebyshev coefficients of U, where U is the representation on the grid. We first compute the chebyshev coefficients rigorously up to order N_ and then compute the norm.

    #### for the second, we use the Lebesgue constant, and take the supremum on the grid.
    S = space(U[1])
    K = order(S)
    npts_= length(U)
    U_cheb = Vector{Sequence{Chebyshev,Vector{Interval{Float64}}}}(undef, 2*K+1)
    for k = -K:K 
         U_cheb[K+k+1] = rifft!(complex.([U[j][k] for j = 1:npts_]), Chebyshev(N_))  #computation of the Chebyshev coefficients
    end

    ΛN = comp_lamb(N_)
    return minimum([norm(norm.(U_cheb,1),1) ΛN*norm(norm.(U,1),Inf)])
end



function norm_seq_grid(U::Vector{Sequence{CartesianProduct{Tuple{ParameterSpace, CosFourier{Interval{Float64}}}}, Vector{Interval{Float64}}}},N_)
    #### computation of an upper bound for the supremum of a seuquence norm. The upper bound is achieved using the minimum of 2 upper bounds.

    ### For the first one, we compute the ℓ1 norm of the chebyshev coefficients of U, where U is the representation on the grid. We first compute the chebyshev coefficients rigorously up to order N_ and then compute the norm.

    #### for the second, we use the Lebesgue constant, and take the supremum on the grid.
    S = space(U[1])
    K = order(S[2])
    npts_= length(U)
    U_cheb = Vector{Sequence{Chebyshev,Vector{Interval{Float64}}}}(undef, K+2)
    for k = 0:K 
         U_cheb[k+2] = rifft!(complex.([component(U[j],2)[k] for j = 1:npts_]), Chebyshev(N_))  #computation of the Chebyshev coefficients
    end
    U_cheb[1] = rifft!(complex.([component(U[j],1)[1] for j = 1:npts_]), Chebyshev(N_))  #computation of the Chebyshev coefficients
    weight = interval.(2*ones(length(U_cheb))) ; weight[1] = interval(1) ; weight[2] = interval(1)

    ΛN = comp_lamb(N_)
    return minimum([norm(norm.(weight.*U_cheb,1),1) ΛN*norm(norm.(U,1),Inf)])
end



function norm_op_grid(A, N_, K1, K2)

    A_cheb = [rifft!(complex.([A[k][i,j] for k = 1:length(A)]), Chebyshev(N_)) for i = indices(codomain(A[1])), j = indices(domain(A[1]))]

    N1, N2 = size(A_cheb)
    A_norm = interval.(zeros(size(A_cheb)[1], size(A_cheb)[2]))

    for i = K1:N1
        for j = K2:N2
            A_norm[i,j] = norm(A_cheb[i,j],1)
        end
    end
    weight1 = interval.(2*ones(N1)) ; weight1[1] = interval(1)
    weight2 = interval.(0.5*ones(N2)) ; weight2[1] = interval(1)
    return opnorm(weight1.*A_norm.*weight2',1)
end




 function cheb2grid(x::VecOrMat{<:Sequence}, N_fft)
          vals = fft.(x, N_fft)
          return [real.(getindex.(vals, i)) for i ∈ eachindex(vals[1])]
      end
      







function proof_branch(K0,N,m0,m1,h,γ,g,U_fft,U_fft_big ; ∂m_control=0)
    d = h 
  N_ = 9N 
N_fft_ = nextpow(2, 2N_ + 1)
npts_ = N_fft_ ÷ 2 + 1


a_fft = [Sequence(CosFourier(K0,interval(1)), component(U_fft[n],2)[:]) for n = 1:2*(npts_-1)]
Q_fft = [component(U_fft[n],1)[1] for n = 1:2*(npts_-1)]

 K1 = 2*K0
 K2 = K1 + 2*K0
  K3 = K1 + 3*K0
 K4 = K1 + 4*K0
# Sc0 = CosFourier(K1,interval(1))⊗Chebyshev(N0)
 S1 = ParameterSpace()×CosFourier(K1,interval(1))
 S2 = ParameterSpace()×CosFourier(K2,interval(1))
 S3 = ParameterSpace()×CosFourier(K3,interval(1))
 S4 = ParameterSpace()×CosFourier(K4,interval(1))

 Sa1 = CosFourier(K1,interval(1))
 Sa2 = CosFourier(K2,interval(1))
 Sa3 = CosFourier(K3,interval(1))
 Sa4 = CosFourier(K4,interval(1))

m_ = [interval(0.5) * (interval(m0) + interval(m1)) + interval(0.5)*cospi(interval(2)*interval(k) / interval(N_fft_)) * (interval(m1) - interval(m0)) for k ∈ 0:npts_-1]
 m_fft_ = [m_ ; reverse(m_)[2:npts_-1]]


# DF_grid0 = DF_grid(m_fft_,h,γ,g,U_fft,CosFourier(K2, interval(1)),2*(nextpow(2, 2N + 1)÷ 2))
# A_grid = inv.([mid.(project(DF_grid0[i],S1,S1)) for i=1:2*(nextpow(2, 2N + 1)÷ 2)])

DF_grid0 = DF_grid(m_fft_,h,γ,g,U_fft,CosFourier(K3, interval(1)),2*npts_-2)

#### we rename npts_ to the number of grid points 2*(npts_-1)
npts_ = 2*(npts_-1)

e0 = Sequence(CosFourier(K0,interval(1)), interval.(zeros(K0+1)))
e0[0] = interval(1)
 
### construction of the v_i fucntions on the grid. First we initialize them
v1_grid = Vector{Sequence{CosFourier{Interval{Float64}},Vector{Interval{Float64}}}}(undef, npts_)
v2_grid = Vector{Sequence{CosFourier{Interval{Float64}},Vector{Interval{Float64}}}}(undef, npts_)
v3_grid = Vector{Sequence{SinFourier{Interval{Float64}},Vector{Interval{Float64}}}}(undef, npts_)
v4_grid = Vector{Sequence{CosFourier{Interval{Float64}},Vector{Interval{Float64}}}}(undef, npts_)

w = Vector{Sequence{Fourier{Interval{Float64}},Vector{Interval{Float64}}}}(undef, npts_)
wc = Vector{Sequence{CosFourier{Interval{Float64}},Vector{Interval{Float64}}}}(undef, npts_)
ws = Vector{Sequence{SinFourier{Interval{Float64}},Vector{Interval{Float64}}}}(undef, npts_)
#### construction of the functions v_i and w. 
for n = 1:npts_
    v1_grid[n] = interval(2)*(m_fft_[n]*Dy_fct_1D(e0,d) + interval(0.5)*γ*Dy_fct_1D((a_fft[n]^2),d) - γ*a_fft[n]*(Dy_fct_1D(a_fft[n],d)))
    v2_grid[n] = interval(2)*g*((Derivative(1)*a_fft[n])^2 + (Dy_fct_1D(a_fft[n],d))^2)
    v3_grid[n] =  -interval(2)*(Q_fft[n] - interval(2)*g*(a_fft[n]))*(Derivative(1)*a_fft[n])
    v4_grid[n] = -interval(2)*(Q_fft[n] - interval(2)*g*(a_fft[n]))*(Dy_fct_1D(a_fft[n],d)) 
    w[n] = mid.(im*sin2exp(v3_grid[n]) + cos2exp(v4_grid[n])) ; w[n] = real.(w[n]) ; w[n] = interval.(mid.(inv(w[n])))
    wc[n] = Sequence(CosFourier(K0,interval(1)),real.((interval.(0.5)*(w[n] + reverse(w[n])))[0:K0]))
    ws[n] = Sequence(SinFourier(K0,interval(1)),real.((interval.(0.5)*(w[n] - reverse(w[n])))[1:K0]))
end

# #### inverse of the derivative (will be used for high frequencies)
# Dx_inv = LinearOperator(SinFourier(K3,interval(1)), Sa4, interval.(zeros(K3+1,K3)))
# for k = 1:K4
#     Dx_inv[k,k] = - interval(1)/interval(k)
# end

##### construction of the approximate inverse  A = A_K - A_K DF A∞ + A∞  (Schur complement, eq (6.34)) #####
### inverse of the derivative on the K2 modes (to build the tail inverse with codomain Sa2)
Dx_inv = LinearOperator(SinFourier(K3,interval(1)), Sa3, interval.(zeros(K3+1,K3)))
for k = 1:K3
    Dx_inv[k,k] = - interval(1)/interval(k)
end
### analytic tail inverse  Â∞ = π^{>K1} D̂y⁻¹ M_w  (eq (6.29)) on the grid, codomain Sa2
A_inf_grid = [Dy_op_inv(Sa3,d)*project(Multiplication(wc[n]),Sa3,Sa3) - Dx_inv*project(Multiplication(ws[n]),Sa3,SinFourier(K3,interval(1))) for n = 1:npts_]
for n = 1:npts_
    A_inf_grid[n][0:K1,:] .= interval.(zeros(K1+1,K3+1)) ### we only keep the tail of the operator
end
### embedding of the tail inverse into H1 :  A∞ U = (0, Â∞ a)  (eq (6.33))
A∞_grid = [ (op = LinearOperator(S3, S3, interval.(zeros(K3+2,K3+2))) ; component(op,2,2) .= coefficients(A_inf_grid[n]) ; op) for n = 1:npts_ ]
### Schur-complement operator  M = DF(U) - DF(U) A∞ DF(U)  on the grid
M_grid = [DF_grid0[n] - DF_grid0[n]*A∞_grid[n]*DF_grid0[n] for n = 1:npts_]

### finite Schur inverse  A_K = π^{≤K1} A_K π^{≤K1}, made rigorous on the grid via Chebyshev interpolation
A_K_grid = inv.([mid.(project(M_grid[i],S1,S1)) for i=1:npts_])
A_cheb = [interval.(mid.(rifft!(complex.([A_K_grid[k][i,j] for k = 1:length(A_K_grid)]), Chebyshev(N)))) for i = indices(codomain(A_K_grid[1])), j = indices(domain(A_K_grid[1]))]
A_K_grid = real.(cheb2grid(A_cheb, N_fft_))
A_K_grid = [LinearOperator(S1,S1 ,A_K_grid[i]) for i = 1:npts_]

######### Computation of the Z1 bound ###########
Id = LinearOperator(S1, S1, interval.(1.0*I[1:K1+2, 1:K1+2]))
### Z10 = residual of  A_K (DF - DF A∞ DF) = A_K M   (eq (6.35))
AKM_grid = A_K_grid.*M_grid
AKM_Z_grid = [Id-D0_op(S1)*project(AKM_grid[i],S1,S1)*D0_inv_op(S1) for i = 1:npts_]
Z10 = norm_op_grid(AKM_Z_grid, N_, 1, 1)
println("Bound Z10 : Z10 = $Z10")

### operator D̂y⁻¹ Mw on the grid (tail, codomain Sa4), used for Z11
Dy_Mw_grid = [Dy_op_inv(Sa4,h)*project(Multiplication(wc[n]),Sa3,Sa4) - Dx_inv*project(Multiplication(ws[n]),Sa3,SinFourier(K4,interval(1))) for n = 1:npts_]
for n = 1:npts_
    Dy_Mw_grid[n][0:K1,:] .= interval.(zeros(K1+1,K3+1)) ### we only keep the tail of the operator
end
### Z11 (unchanged) : Dy^{-1} Mw DF on the tail rows
CDF_grid = [D0_op(Sa4)*Dy_Mw_grid[n]*project(component(DF_grid0[n],2,2),Sa1,Sa3)*D0_inv_op(Sa1) for n = 1:npts_]
Z11 = norm_op_grid(CDF_grid, N_, K1+2, 1)
println("Bound Z11 : Z11 = $Z11")

### Z12 = ‖A_K DF(U) π^{>K1}‖   (eq (6.35), plays the role of the former Z14)
AKDF_grid = A_K_grid.*DF_grid0
AKDF_Z_grid = [D0_op(S1)*AKDF_grid[i]*D0_inv_op(S3) for i = 1:npts_]
Z12 = norm_op_grid(AKDF_Z_grid, N_, 1, K1+2)
println("Bound Z12 : Z12 = $Z12")


##### Computation of the different components of Z∞
### We start by computing ϵ
w_prod = [-ws[n]*v3_grid[n] + wc[n]*v4_grid[n] - interval(1) for n = 1:npts_]
w_v4 = [w[n]*cos2exp(v4_grid[n]) for n=1:npts_]

ϵ = norm_seq_grid(w_prod, N_) +  interval(2)*norm_seq_grid(w_v4, N_)*exp(-interval(2*K0)*h)/(interval(1) - exp(-interval(2*K0)*h))
ϵ = ϵ*interval(K1+2)/interval(K1+1) 


D = LinearOperator(Fourier(K0,interval(1)),Fourier(K0,interval(1)), interval.(1.0*I[1:2K0+1, 1:2K0+1]))
for k = -K0:-1
    D[k,k] = interval(k) - abs(interval(k))*coth(d*abs(interval(k)))
end
for k = 1:K0
    D[k,k] = interval(k) - abs(interval(k))*coth(d*abs(interval(k)))
end
    D[0,0] = - interval(1)/d

#### first term of Z∞    
w_prod = [w[n]*(cos2exp(v1_grid[n])*(D*cos2exp(a_fft[n])) + cos2exp(v2_grid[n])) for n=1:npts_]
Z∞ = norm_seq_grid(w_prod,N_)/(interval(K1+1))^2*interval(K1+2)
### adding the second term
w_prod = [w[n]*cos2exp(v1_grid[n]) for n=1:npts_]
a_fft_D = [D0_op(Sa1)*a_fft[n] for n = 1:npts_]

Z∞ = Z∞ + interval(4)*norm_seq_grid(w_prod,N_)*exp(-interval(2*K0)*d)/(interval(1) - exp(-interval(2*K0)*d))*interval(K0)*interval(K1+2)/interval(K1+1)*norm_seq_grid(a_fft_D,N_)
Z∞ = Z∞ + ϵ #### we finish by adding epsilon computed above

println("Bound Z∞ : Z∞ = $Z∞")

Z1 = maximum([Z10+Z11 (interval(1)+Z12)*Z∞])
println("Bound Z1 : Z1 = $Z1")

########## Computation of the Y bound ###########
#### When the candidate U_fft_big is given in extended precision (Interval{BigFloat}), the residual
#### F(U,m) is evaluated in BigFloat (so the cancellation is resolved) and then rigorously enclosed
#### in Interval{Float64}. All the remaining operations stay in standard precision.
if eltype(coefficients(U_fft_big[1])) == Interval{BigFloat}
    np = N_fft_ ÷ 2 + 1
    m_big = [interval(big(0.5))*(interval(big(m0)) + interval(big(m1))) + interval(big(0.5))*cospi(interval(big(2))*interval(big(j))/interval(big(N_fft_)))*(interval(big(m1)) - interval(big(m0))) for j ∈ 0:np-1]
    m_fft_big = [m_big ; reverse(m_big)[2:np-1]]
    h_big = interval(big(inf(h)),big(sup(h))) ; γ_big = interval(big(inf(γ)),big(sup(γ))) ; g_big = interval(big(inf(g)),big(sup(g)))
    F_grid0 = F_grid(m_fft_big, h_big, γ_big, g_big, U_fft_big, S4, npts_)
    F_grid0 = [interval.(Float64.(inf.(F_grid0[i]),RoundDown), Float64.(sup.(F_grid0[i]),RoundUp)) for i = 1:npts_]
else
    F_grid0 = F_grid(m_fft_, h, γ, g, U_fft_big, S4, npts_)
end
### finite component :  A_K (I - DF(U) A∞) F(U)
AKF_grid = [A_K_grid[i]*(F_grid0[i] - DF_grid0[i]*(A∞_grid[i]*F_grid0[i])) for i = 1:npts_]
Y = norm_seq_grid([D0_fct(AKF_grid[i]) for i = 1:npts_],N_)
### tail component :  (K1+2)/(K1+1) ( ‖π^{>K1}(w*f)‖ + ‖π^{<-K1}(w⋆*f)‖ )
wf  = [ (s = w[n]*cos2exp(component(F_grid0[n],2)) ; Nf = order(s) ; c = copy(coefficients(s)) ; c[1:Nf+K1+1] .= interval(0) ; Sequence(space(s), c)) for n = 1:npts_ ]
wfc = [ (s = reverse(w[n])*cos2exp(component(F_grid0[n],2)) ; Nf = order(s) ; c = copy(coefficients(s)) ; c[Nf-K1+1:2*Nf+1] .= interval(0) ; Sequence(space(s), c)) for n = 1:npts_ ]
Y = Y + interval(K1+2)/interval(K1+1)*(norm_seq_grid(wf,N_) + norm_seq_grid(wfc,N_))
println("Bound Y : Y = $Y")

######## computation Z2 bound ##############

r0 = interval(1e-4)


norm_ell_1_a = norm_seq_grid(a_fft,N_)
norm_H1_a = norm_seq_grid(a_fft_D,N_)

####### computation norm of A #########
### full approximate inverse  A = A_K - A_K DF(U) A∞ + A∞  on the grid  (S1 -> S2)
A_grid = [project(A_K_grid[i],S1,S1) - project(A_K_grid[i]*DF_grid0[i]*A∞_grid[i],S1,S3) + project(A∞_grid[i],S1,S2) for i = 1:npts_]
norm_A = norm_op_grid([D0_op(S2)*A_grid[i] for i = 1:npts_], N_, 1, 1)
norm_w = norm_seq_grid(w,N_)
norm_A = maximum([norm_A  (interval(1)+Z12)*interval(K1+2)/interval(K1+1)*norm_w])

######## computation of ZD2 ##############

fct1 = [m_fft_[n]/d + interval(0.5)*γ*Dy_fct_1D((a_fft[n]^2),d) - γ*a_fft[n]*(Dy_fct_1D(a_fft[n],d)) for n = 1:npts_]
fct2 = [Q_fft[n]-interval(2)*g*(a_fft[n]) for n = 1:npts_]

ZD2 = interval(2)*γ^2*(interval(2)*norm_H1_a + norm_ell_1_a)^2/d^2 
ZD2 += interval(6)*abs(γ)/d*norm_seq_grid(fct1,N_)
ZD2 += (interval(1) + interval(1)/d^2)*(interval(4)*(interval(1)+interval(2)*g)*norm_H1_a + norm_seq_grid(fct2,N_))

###### computation of ZD3 #############
max_Q = norm(rifft!(complex.(Q_fft), Chebyshev(N_)),1) 

ZD3 = r0*(interval(27)*γ^2/d^2*(norm_H1_a+r0) + interval(4)*g*(interval(1) + interval(1)/d^2) )
ZD3 += r0*(abs(max_Q) +interval(2)*g)*(interval(1) + interval(1)/d^2) 
ZD3 += r0*interval(3)*(abs(max_Q)+r0)*(interval(1) + interval(1)/d^2)

Z2 = norm_A*(ZD2 + ZD3)
println("Bound Z2 : Z2 = $Z2")





if inf((interval(1)-Z1)^2) > sup(interval(2)*Z2*Y)
    rm = (interval(1)-interval(sup(Z1)) - sqrt((interval(1)-interval(sup(Z1)))^2-interval(2)*Z2*Y))/Z2
    rp = (interval(1)-interval(sup(Z1)) + sqrt((interval(1)-interval(sup(Z1)))^2-interval(2)*Z2*Y))/Z2
    rp = minimum([rp (interval(1)-Z1)/Z2 r0])
    if (sup(abs(rm))<inf(r0))&&(sup(Z1 + Z2*rm) < 1)
        ## if ∂m_control=0, we simply return the radiae of the proof and the bounds. Otherwise, we return the enclosure of ∂mŨ(m).
        if ∂m_control == 0
            println("Proof successful with radius in : [$rm , $rp ]")
            return (rm=rm, rp=rp, Z1=Z1, Z2=Z2, normA=norm_A, Y=Y)
        else 

            ##In order to prove that the branch undergoes a transition from graph to overhanging, we need to control ∂mŨ(m). 
            ## If the proof succeeds, we will return an approximation V(m) for ∂mŨ(m) as well as a rigorous error bound on the difference \|V(m) - ∂mŨ(m)\|_X. The error bound is given in Section 9

            ### we first compute V(m) = ∂m Ū(m) on the grid.
            U_cheb = Vector{Sequence{Chebyshev,Vector{Interval{Float64}}}}(undef, K0+2)
            for n = 1:K0+2
                U_cheb[n] = real.(rifft!(complex.([coefficients(U_fft[j])[n] for j = 1:npts_]), Chebyshev(N)))  #computation of the Chebyshev coefficients
            end
            ## the Chebyshev coefficients are in τ ∈ [-1,1] and m is affine in τ, so ∂m = (dτ/dm) ∂τ , dτ/dm = 2/(m1-m0)
            dτdm = interval(2)/(interval(m1) - interval(m0))
            Vm_cheb = [dτdm*(Derivative(1)*U_cheb[n]) for n = 1:K0+2]
            Vm_vec  = cheb2grid(Vm_cheb, N_fft_)
            ## wrap every grid node as a sequence in ParameterSpace × CosFourier(K0)
            Vm_grid = [Sequence(ParameterSpace()×CosFourier(K0,interval(1)), Vm_vec[n]) for n = 1:npts_]
            ## same data embedded into CosFourier(K3) (the domain of DF_grid0) for the residual below
            Vm_K3 = [ (s = Sequence(ParameterSpace()×CosFourier(K3,interval(1)), interval.(zeros(K3+2))) ; coefficients(s)[1:K0+2] .= Vm_vec[n] ; s) for n = 1:npts_ ]
        
            ## computation of the defect 
            ## D_mF has no Q-component and equals (2/h)·v1 in the a-component ; build it as a
            ## ParameterSpace × CosFourier(K3) element so the sum with DF_grid0·V is well defined.
            DmF = [ (s = Sequence(ParameterSpace()×CosFourier(K3,interval(1)), interval.(zeros(K3+2))) ;
                     component(s,2) .= coefficients(project(I"2"/h*(m_fft_[n]/h + γ*I"0.5"*Dy_fct_1D(a_fft[n]^2,h) - γ*a_fft[n]*Dy_fct_1D(a_fft[n],h)), CosFourier(K3,interval(1)))) ; s) for n = 1:npts_ ]
            err_Vm = [DF_grid0[n]*Vm_K3[n] + DmF[n] for n = 1:npts_]

            ### computation of the first component of ϵ, which is \| A(D_UF*V(m) + D_m F) \|_X. We compute it exactly as the Y bound
            A_err_Vm = [A_K_grid[i]*(err_Vm[i] - DF_grid0[i]*(A∞_grid[i]*err_Vm[i])) for i = 1:npts_]
            ϵ = norm_seq_grid([D0_fct(A_err_Vm[i]) for i = 1:npts_],N_)
            wf  = [ (s = w[n]*cos2exp(component(err_Vm[n],2)) ; Nf = order(s) ; c = copy(coefficients(s)) ; c[1:Nf+K1+1] .= interval(0) ; Sequence(space(s), c)) for n = 1:npts_ ]
            wfc = [ (s = reverse(w[n])*cos2exp(component(err_Vm[n],2)) ; Nf = order(s) ; c = copy(coefficients(s)) ; c[Nf-K1+1:2*Nf+1] .= interval(0) ; Sequence(space(s), c)) for n = 1:npts_ ]
            ϵ = ϵ + interval(K1+2)/interval(K1+1)*(norm_seq_grid(wf,N_) + norm_seq_grid(wfc,N_))

            ### the isolation radius is rm =  r⁻, NOT the a-priori r0 = 1e-4 used for Z2.
            ### ϵ = [ ‖A·residual‖_X  +  Z2 rm ‖V‖_X  +  (parameter-direction Lipschitz) ] / (1 - Z1 - Z2 rm)
            ϵ = ϵ + Z2*rm*norm_seq_grid([D0_fct(Vm_grid[i]) for i = 1:npts_],N_)        ### Lipschitz term, U-direction
            ϵ = ϵ + norm_A*I"4"/h*abs(γ)*coth(h)*(norm_H1_a + rm)*rm                    ### Lipschitz term, parameter direction
            ϵ = ϵ/(I"1" - Z1 - Z2*rm)                                                   ### Neumann series factor (denominator)

            return (rm=rm, rp=rp, Z1=Z1, Z2=Z2, normA=norm_A, Y=Y), Vm_grid, ϵ
        end
    else
        display("second condition not verified")
        return naN
    end
else
    display("first condition not verified")
    return naN
end  

end







#### Extended-precision computation of the branch. For each segment [mi0, mi] the Newton method
#### is run in BigFloat (≈ `digits` decimal digits) down to residual ≤ `tol`, and the grid
#### representation U_fft_grid is saved (as Interval{BigFloat}) to a JLD2 file named by the m-interval.
#### No proof is run here: load the files afterwards to run proof_branch.
function save_branch(m0, m1, dm, N, U0, K0, γ, g, h, k; digits=128, tol=1e-30, folder=".")

    setprecision(digits)   #### ≈ `digits` decimal digits of precision

    N_ = 9N
    N_fft_ = nextpow(2, 2N_ + 1)
    npts_ = N_fft_ ÷ 2 + 1

    Tb = Sequence{CartesianProduct{Tuple{ParameterSpace, CosFourier{Int64}}}, Vector{BigFloat}}

    #### parameters and initial guess in extended precision
    γb = big(mid(γ)) ; gb = big(mid(g)) ; hb = big(mid(h)) ; kb = big(mid(k))
    U_guess = Sequence(ParameterSpace()×CosFourier(K0,1), big.(coefficients(mid.(U0))))

    for mi = m0-dm : -dm : m1
        mi0 = mi + dm
        println("Extended-precision solve on segment  m ∈ [$mi , $mi0]")

        #### Chebyshev nodes in m (extended precision)
        m = [big(0.5)*(big(mi0) + big(mi)) + big(0.5)*cospi(big(2)*big(j)/big(N_fft_))*(big(mi) - big(mi0)) for j ∈ 0:npts_-1]

        #### Newton at each node, warm-started node by node
        Uc_fft = Vector{Tb}(undef, npts_)
        Uc_fft[npts_] = Newton_method(m[npts_], γb, kb, gb, hb, U_guess, 100, tol)
        for n = npts_-1:-1:1
            Uc_fft[n] = Newton_method(m[n], γb, kb, gb, hb, Uc_fft[n+1], 100, tol)
        end
        res = maximum([Float64(norm(F(m[n], γb, kb, gb, hb, Uc_fft[n]),1)) for n = 1:npts_])
        println("   max Newton residual on the segment ≈ $res")

        #### Chebyshev representation and grid representation (extended precision)
        U_fft = [Uc_fft ; reverse(Uc_fft)[2:npts_-1]]
        U_cheb = Vector{Sequence{Chebyshev,Vector{BigFloat}}}(undef, K0+2)
        ### we compute the Chebyshev coefficients of a polynomial of order N
        for n = 1:K0+2
            U_cheb[n] = real.(rifft!(complex.([coefficients(U_fft[j])[n] for j = 1:2*(npts_-1)]), Chebyshev(N)))  #computation of the Chebyshev coefficients
        end
        ## evaluate rigorously our polynomial of order N on the 9N grid
        U_fft_grid = cheb2grid(U_cheb, N_fft_)
        U_fft_grid = [Sequence(ParameterSpace()×CosFourier(K0,interval(1)), interval.(U_fft_grid[i])) for i = 1:2*(npts_-1)]

        #### save the segment ; the m-interval is in the file name and also stored in the file
        fname = joinpath(folder, "U_fft_grid_m_$(round(mi0,digits=6))_$(round(mi,digits=6)).jld2")
        save(fname, "U_fft_grid", U_fft_grid, "mi0", mi0, "mi", mi, "K0", K0, "N", N)
        println("   saved  →  $fname")

        U_guess = Uc_fft[1]   #### warm start for the next segment (shared endpoint m = mi)
    end
end


