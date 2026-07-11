function cos2exp(a)
    N = size(a)[1]-1
    return Sequence(Fourier(N,interval(1)), vec([reverse(a[1:N]); a[0:N]]))
end

function sin2exp(a)
    N = size(a)[1]
    return Sequence(Fourier(N,interval(1)), vec([im*reverse(a[1:N]); interval(0) ; -im*a[1:N]]))
end


function Dy_op_inv(Sc,h::Interval{Float64})
    K0 = order(Sc)
    dn = interval.(1:K0)
    dn = [h ; interval(1) ./ dn]
    return  LinearOperator(Sc, Sc, Diagonal( dn ))
end 





function Dy_op_1D(Sc,h::Interval{Float64})
    K = order(Sc)
    dn = interval.(1:K)
    dn = [interval(1)/h ; dn.*coth.(h*dn)] 

    return LinearOperator(Sc, Sc, Diagonal( dn ))
end 


function Dy_op_1D(Sc,h::Float64)
    K = order(Sc)
    dn = (1:K)
    dn = [(1)/h ; dn.*coth.(h*dn)] 

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




function Dy_fct_1D(U,h::Interval{Float64})
    Sc = space(U)
    K = order(Sc)
    dn = interval.(1:K)
    dn = [interval(1)/h ; dn.*coth.(h*dn)]
    return dn.*U
end 


function Dy_fct_1D(U,h::Interval{BigFloat})
    Sc = space(U)
    K = order(Sc)
    dn = interval.(big.(1:K))
    dn = [interval(big(1))/h ; dn.*coth.(h*dn)]
    return dn.*U
end 

function Dy_fct_1D(U,h::BigFloat)
    Sc = space(U)
    K = order(Sc)
    dn = big.(1:K)
    dn = [big(1)/h ; dn.*coth.(h*dn)]
    return dn.*U
end 



function Dy_fct_1D(U,h::Float64)
    Sc = space(U)
    K = order(Sc)
    dn = (1:K)
    dn = [(1)/h ; dn.*coth.(h*dn)]
    return dn.*U
end 





function G_grid(ϵ,h,γ,m̄,g,U,S,npts_)
    T = CartesianProduct{Tuple{ParameterSpace, CosFourier{Interval{Float64}}}}
    F_grid = Vector{Sequence{T,Vector{Interval{Float64}}}}(undef, npts_)
    F_raw  = Vector{Sequence{T,Vector{Interval{Float64}}}}(undef, npts_)   #### NEW (Lemma 7.4): raw residual f = F(V,m̄−ϵ²)
    a = [component(U[n],2) for n = 1:npts_]
   
    Q = [component(U[n],1)[1] for n = 1:npts_]
    e0 = Sequence(space(a[1]), interval.(zeros(size(a[1])[1]))) 
    e0[0] = interval(1.0)
    m = m̄ .- ϵ.^2 

    n0 = Int(npts_ /2) + 1
    J = [1:n0-1 ; n0+1:npts_]

    for n = J
        F_grid0 = Sequence(S, interval.(zeros(order(S[2])+2)))
        component(F_grid0,1)[1] .= a[n][0] - h
        component(F_grid0,2) .= project( (m[n]*Dy_fct_1D(e0,h) + interval(0.5)*γ*Dy_fct_1D(a[n]^2,h) - γ*a[n]*(Dy_fct_1D(a[n],h)))^2 - (Q[n] - interval(2)*g*a[n])*((Derivative(1)*a[n])^2 + Dy_fct_1D(a[n],h)^2), S[2])
        F_raw[n] = F_grid0                                              #### NEW (Lemma 7.4): store raw residual (component 2 is f)
        F_grid[n] = interval(1)/ϵ[n]^2*op_P(ϵ[n],ParameterSpace()×S[2])*F_grid0
    end
    
    F_grid[n0] = Ḡ(project(U[n0],S), S, m̄ , γ, g, h)

    F_raw[n0]  = Sequence(S, interval.(zeros(order(S[2])+2)))           #### NEW (Lemma 7.4): flat state is an exact solution, so f = 0 at the bifurcation

    return F_grid, F_raw
end



function DG_grid(ϵ,h,γ,m̄,g,U,S,npts_)
    T = CartesianProduct{Tuple{ParameterSpace, CosFourier{Interval{Float64}}}}
    DF_grid = Vector{LinearOperator{T,T,Matrix{Interval{Float64}}}}(undef, npts_)
    a = [component(U[n],2) for n = 1:npts_]
  
    Q = [component(U[n],1)[1] for n = 1:npts_]
    m = m̄ .- ϵ.^2

    n0 = Int(npts_ /2) + 1
    J = [1:n0-1 ; n0+1:npts_]

    for n ∈ J
        Ma = project(Multiplication(a[n]), S, S)
        Dy = Dy_op_1D(S,h)
        C = γ*(Dy*Ma - Ma*Dy - project(Multiplication(Dy_fct_1D(a[n],h)), S, S))

        K = order(S)
        Dx = project(Derivative(1), S, SinFourier(K,interval(1)))

        e0 = Sequence(space(a[n]), interval.(zeros(size(a[n])[1])))
        e0[0] = interval(1)

        v1 = interval(2)*(m[n]*Dy*e0 + interval(0.5)γ*Dy_fct_1D(a[n]^2,h) - γ*a[n]*Dy_fct_1D(a[n],h))
        v2 = interval(2)*g*((Derivative(1)*a[n])^2 + (Dy_fct_1D(a[n],h))^2)
        v3 = -interval(2)*(Q[n] - interval(2)*g*(a[n]))*(Derivative(1)*a[n])
        v4 = -interval(2)*(Q[n] - interval(2)*g*(a[n]))*(Dy_fct_1D(a[n],h)) 

        DF0 = LinearOperator(ParameterSpace()×S, ParameterSpace()×S, interval.(zeros(K+2, K+2)))
        e0 = interval.(zeros(1,K+1)) ; e0[1] = interval(1.0)
        component(DF0,1,2).= e0
        component(DF0,2,1) .= - project((Dx*a[n])^2 + (Dy_fct_1D(a[n],h))^2, S)
        component(DF0,2,2) .= project(Multiplication(v1), S, S)*C + project(Multiplication(v2), S, S) + project(Multiplication(v3), SinFourier(K,1), S)*Dx + project(Multiplication(v4), S, S)*Dy
        DF_grid[n] = op_P(ϵ[n],ParameterSpace()×S)*DF0*op_P(ϵ[n],ParameterSpace()×S)
    end

    #### W at the bifurcation point
    S0 = ParameterSpace()×S
    DF_grid[n0] = DḠ(project(U[n0],S0), S0 , m̄ , γ, g, h)
    return DF_grid
end




######### zero finding problem construction 

function op_P(ϵ::Interval{Float64},S)
    K1 = order(S[2])
    P = LinearOperator(S, S, interval.(1.0*I[1:K1+2, 1:K1+2]))
    component(P,2,2)[1,1] .= interval(1)/ϵ
    return P
end


function op_P(ϵ::Interval{BigFloat},S)
    K1 = order(S[2])
    P = LinearOperator(S, S, interval.(big.(1.0*I[1:K1+2, 1:K1+2])))
    component(P,2,2)[1,1] .= interval(big(1.0))/ϵ
    return P
end

function op_P(ϵ::Float64,S)
    K1 = order(S[2])
    P = LinearOperator(S, S, 1.0*I[1:K1+2, 1:K1+2])
    component(P,2,2)[1,1] .= 1/ϵ
    return P
end


function F(m,γ,g,h,U)

    Q = component(U,1)[1]
    a = component(U,2)

    e0 = Sequence(space(a), zeros(size(a)[1]))
    e0[0] = 1.0
    F0 = Sequence(space(U), zeros(size(U)[1]))

    component(F0,1)[1] .= a[0] - h
    component(F0,2) .= project( (m*Dy_fct_1D(e0,h) + γ/2*Dy_fct_1D(a^2,h) - γ*a*(Dy_fct_1D(a,h)))^2 - (Q-2*g*a)*((Derivative(1)*a)^2 + (Dy_fct_1D(a,h))^2) , space(a))
    return F0
end



function DF(m::Float64,γ,g,h,U)

      Q = component(U,1)[1]
    a = component(U,2)
    N = size(a)[1]-1 ; Sa = space(a)
    n = 1:N 
Dy = LinearOperator(Sa, Sa, Diagonal( vec([1/h ; n.*coth.(h*n)]) ))
Dx = project(Derivative(1), Sa, SinFourier(N,1))

    N = order(space(a))
    Ma = project(Multiplication(a), space(a), space(a))
    C = γ*(Dy*Ma - Ma*Dy - project(Multiplication(Dy_fct_1D(a,h)), space(a), space(a)))

    e0 = Sequence(space(a), zeros(size(a)[1]))
    e0[0] = 1.0

    v1 = 2*(m*Dy*e0 + γ/2*Dy*(a^2) - γ*a*(Dy*a))
    v2 = 2*g*((Dx*a)^2 + (Dy*a)^2)
    v3 = -2*(Q - 2*g*(a))*(Dx*a)
    v4 = -2*(Q - 2*g*(a))*(Dy*a) 

    DF = LinearOperator(space(U), space(U), zeros(size(U)[1], size(U)[1]))
    e0 = zeros(1,N+1) ; e0[1] = 1.0
    component(DF,1,2).= e0
    component(DF,2,2).= project(Multiplication(v1), space(a), space(a))*C + project(Multiplication(v2), space(a), space(a)) + project(Multiplication(v3), SinFourier(N,1), space(a))*Dx + project(Multiplication(v4), space(a), space(a))*Dy
    component(DF,2,1) .= - project((Dx*a)^2 + (Dy*a)^2, space(a)) 
    return DF
end



function DF(m,γ,g,h,S,U)

    Sa = S[2] ; N = order(Sa)

    Q = component(U,1)[1]
    a = project(component(U,2),Sa)

    n = interval.(1:N) 
    Dy = LinearOperator(Sa, Sa, Diagonal( vec([interval(1)/h ; n.*coth.(h*n)]) ))
    Dx = project(Derivative(1), Sa, SinFourier(N,interval(1)))

    Ma = project(Multiplication(a), Sa, Sa)
    C = γ*(Dy*Ma - Ma*Dy - project(Multiplication(Dy_fct_1D(a,h)), Sa, Sa))

    e0 = Sequence(Sa, interval.(zeros(N+1)))
    e0[0] = interval(1.0)

    v1 = interval(2)*(m*Dy*e0 + γ*interval(0.5)*Dy*(a^2) - γ*a*(Dy*a))
    v2 = interval(2)*g*((Dx*a)^2 + (Dy*a)^2)
    v3 = -interval(2)*(Q - interval(2)*g*(a))*(Dx*a)
    v4 = -interval(2)*(Q - interval(2)*g*(a))*(Dy*a) 

    DF = LinearOperator(S, S, interval.(zeros(N+2, N+2)))
    
    component(DF,1,2) .= coefficients(e0)'
    component(DF,2,2).= project(Multiplication(v1), Sa, Sa)*C + project(Multiplication(v2), Sa, Sa) + project(Multiplication(v3), SinFourier(N,1), Sa)*Dx + project(Multiplication(v4), Sa, Sa)*Dy
    component(DF,2,1) .= - project((Dx*a)^2 + (Dy*a)^2, Sa) 
    return DF
end




function G(U, ϵ, m̄ , γ, g, h)

    a = component(U,2) ; Sa = space(a) ; K0 = order(Sa)

    m = m̄-ϵ^2
    Q̄ = 2*g*h + (m/h - γ*h*0.5)^2
    #### flat solution at the bifurcation point
    ā = Sequence(Sa, zeros(K0+1)) ; ā[0] = h
    Ū = Sequence(space(U),[Q̄ ; ā[0:K0]])

    Pϵ = op_P(ϵ, space(U))
    if ϵ == 0
        return mid.(Ḡ(Ū, space(U), m̄ , γ, g, h))
    else 
        return 1/ϵ^2*Pϵ*F(m,γ,g,h,Ū + ϵ^2*Pϵ*U)
    end
end



function DG(U, ϵ, m̄ , γ, g, h)

    a = component(U,2) ; Sa = space(a) ; K0 = order(Sa)

    m = m̄-ϵ^2
    Q̄ = 2*g*h + (m/h - γ*h*0.5)^2
    #### flat solution at the bifurcation point
    ā = Sequence(Sa, zeros(K0+1)) ; ā[0] = h
    S = ParameterSpace()×CosFourier(K0,1) 
    Ū = Sequence(S,[mid(Q̄) ; mid.(ā[0:K0])])

    Pϵ = op_P(ϵ, space(U))
    if ϵ == 0 
        return mid.(DḠ(Ū, space(U), m̄ , γ, g, h))
    else 
        return Pϵ*DF(m,γ,g,h,Ū + ϵ^2*Pϵ*U)*Pϵ
    end
end





function D2F(m,γ,g,h,S,U,V,W)

    D2F0 = Sequence(S, interval.(zeros(order(S[2])+2)))
    Q = component(U,1)[1]
    a = component(U,2)

    β1 = component(V,1)[1]
    v = component(V,2)

    β2 = component(W,1)[1]
    w = component(W,2)

    K = order(S[2])
    Dx = project(Derivative(1), CosFourier(K,interval(1)), SinFourier(K,interval(1)))

    e0 = interval.(zeros(K+1)) ; e0[1] = interval(1.0) ; e0 = Sequence(S[2], e0)

    component(D2F0,2) .= project( interval(2)*γ^2*(Dy_fct_1D(a*v,h) - a*Dy_fct_1D(v,h) - v*Dy_fct_1D(a,h))*(Dy_fct_1D(a*w,h) - a*Dy_fct_1D(w,h) - w*Dy_fct_1D(a,h)) + interval(2)*γ*(m*Dy_fct_1D(e0,h) + interval(0.5)γ*Dy_fct_1D(a^2,h) - γ*a*(Dy_fct_1D(a,h)))*(Dy_fct_1D(v*w,h) - w*Dy_fct_1D(v,h) - v*Dy_fct_1D(w,h))  + (interval(4)*g*v-interval(2)*β1)*((Dx*a)*(Dx*w)   + (Dy_fct_1D(a,h))*(Dy_fct_1D(w,h))) + (interval(4)*g*w-interval(2)*β2)*((Dx*a)*(Dx*v)   + (Dy_fct_1D(a,h))*(Dy_fct_1D(v,h))) - interval(2)*(Q - interval(2)*g*(a))*((Dx*v)*(Dx*w) + Dy_fct_1D(v,h)*Dy_fct_1D(w,h)), S[2] )

    return D2F0
end




function D2F_op(m,γ,g,h,S,U,V)

    Sa = S[2] ; K = order(Sa)
    D2F0 = LinearOperator(S, S, interval.(zeros(K+2, K+2)))

    Q = component(U,1)[1]
    a = component(U,2)

    β = component(V,1)[1]
    v = component(V,2)

    Sa = S[2] ; K = order(Sa)

    n = interval.(1:K) 
    Dy = LinearOperator(Sa, Sa, Diagonal( vec([interval(1)/h ; n.*coth.(h*n)]) ))
    Dx = project(Derivative(1), Sa, SinFourier(K,interval(1)))

    e0 = interval.(zeros(K+1)) ; e0[1] = interval(1.0) ; e0 = Sequence(S[2], e0)

    Ma = project(Multiplication(a), Sa, Sa)
    MDya = project(Multiplication(Dy_fct_1D(a,h)), Sa, Sa)
    MDxa = project(Multiplication(Dx*a), SinFourier(K,interval(1)), Sa )
    MQa = project(Multiplication(Q - interval(2)*g*a), Sa, Sa)
    MCa = project(Multiplication(Dy_fct_1D(a*v,h) - a*Dy_fct_1D(v,h) - v*Dy_fct_1D(a,h)), Sa, Sa)
    MC2a = project(Multiplication(m*Dy_fct_1D(e0,h) + interval(0.5)γ*Dy_fct_1D(a^2,h) - γ*a*(Dy_fct_1D(a,h))), Sa, Sa)

    Mv = project(Multiplication(v), Sa, Sa)
    MDyv = project(Multiplication(Dy_fct_1D(v,h)), Sa, Sa)
    MDxv = project(Multiplication(Dx*v), SinFourier(K,interval(1)), Sa)
    MQv = project(Multiplication(β - interval(2)*g*v), Sa, Sa)
    
    MDaDv = project(Multiplication((Dx*a)*(Dx*v)   + (Dy_fct_1D(a,h))*(Dy_fct_1D(v,h))), Sa, Sa)

   

    component(D2F0,2,2) .= project( interval(2)*γ^2*MCa*(Dy*Ma - Ma*Dy - MDya) + interval(2)*γ*MC2a*(Dy*Mv - MDyv - Mv*Dy) + interval(4)*g*MDaDv - interval(2)*MQv*(MDya*Dy + MDxa*Dx)  - interval(2)*MQa*(MDyv*Dy + MDxv*Dx) , Sa, Sa)

    component(D2F0,2,1) .= project(-interval(2)*((Dx*a)*(Dx*v)   + (Dy_fct_1D(a,h))*(Dy_fct_1D(v,h))) ,Sa)

    return D2F0
end



function DmF(m,γ,h,U)

    S = space(U)
    DmF0 = Sequence(S, zeros(size(U)[1]))
    a = component(U,2)

    e0 = zeros(size(a)[1]) ; e0[1] = 1.0 ; e0 = Sequence(space(a), e0)

    component(DmF0,2) .= project( 2*(m*Dy_fct_1D(e0,h) + γ/2*Dy_fct_1D(a^2,h) - γ*a*(Dy_fct_1D(a,h)))*Dy_fct_1D(e0,h) , space(a) )

    return DmF0
end

function DmDUF(γ,h,S,U,V)

    DmF0 = Sequence(S, interval.(zeros(order(S[2])+2)))
    a = component(U,2)
    v = component(V,2)

    e0 = interval.(zeros(size(a)[1])) ; e0[1] = interval(1.0) ; e0 = Sequence(space(a), e0)

    component(DmF0,2) .= project( interval(2)*γ*Dy_fct_1D(e0,h)*(Dy_fct_1D(a*v,h) - a*(Dy_fct_1D(v,h)) - v*(Dy_fct_1D(a,h))), S[2] )
    return DmF0
end



function D3UF(γ,h,S,U,V)

    DmF = Sequence(S, interval.(zeros(order(S[2])+2)))
    a = component(U,2)

     β = component(V,1)[1]
    v = component(V,2)

    e0 = interval.(zeros(size(a)[1])) ; e0[1] = interval(1.0) ; e0 = Sequence(space(a), e0)

    component(DmF,2) .= project( interval(6)*γ^2*(Dy_fct_1D(a*v,h) - a*(Dy_fct_1D(v,h)) - v*(Dy_fct_1D(a,h)))*(Dy_fct_1D(v^2,h) - interval(2)*v*Dy_fct_1D(v,h)) - interval(6)*(β - interval(2)*g*v)*(Dy_fct_1D(v,h)^2 + (Derivative(1)*v)^2), S[2])

    return DmF
end



function Ḡ(U, S, m̄ , γ, g, h)
    
    a = component(U,2) 
    Sa = S[2] ; K = order(Sa)
    a = project(a, Sa)
   
    Z̄ = Sequence(S,interval.(zeros(K+2))) ; component(Z̄,2)[1] .= interval(1)
    α = Z̄ .* U ; W = U - α
    α = coefficients(α)

    Q̄ = interval(2)*g*h + (m̄/h - γ*h*interval(0.5))^2
    #### flat solution at the bifurcation point
    ā = Sequence(Sa, interval.(zeros(K+1))) ; ā[0] = h
    Ū = Sequence(S, [Q̄ ; coefficients(ā)])

    DmŪ = Sequence(S, [interval(2)/h*(m̄/h - γ*h*interval(0.5)) ; interval.(zeros(K+1))])
    
    G0 = α .* D2F(m̄,γ,g,h,S,Ū,Z̄,W-DmŪ) - α .* DmDUF(γ,h,S,Ū,Z̄) + interval(1)/interval(6)*(α.^3) .* D3UF(γ,h,S,Ū,Z̄) + DF(m̄,γ,g,h,S,Ū)*W + a[1]^2*interval(0.5)*D2F(m̄,γ,g,h,S,Ū,Z̄,Z̄)

    return project(G0,S)
end




function DḠ(U, S, m̄ , γ, g, h)
    
    a = component(U,2) 
    Sa = S[2] ; K = order(Sa)
    a = project(a, Sa)
   
    Z̄ = Sequence(S,interval.(zeros(K+2))) ; component(Z̄,2)[1] .= interval(1)
    α = Z̄ .* U ; W = U - α
    α = coefficients(α)


    Q̄ = interval(2)*g*h + (m̄/h - γ*h*interval(0.5))^2
    #### flat solution at the bifurcation point
    ā = Sequence(Sa, interval.(zeros(K+1))) ; ā[0] = h
    Ū = Sequence(S, [Q̄ ; coefficients(ā)])

    DmŪ = Sequence(S, [interval(2)/h*(m̄/h - γ*h*interval(0.5)) ; interval.(zeros(K+1))])

    DG0 =  LinearOperator(S, S, interval.(zeros(K+2, K+2)))

    DαG = Z̄ .* D2F(m̄,γ,g,h,S,Ū,Z̄,W-DmŪ) - Z̄ .* DmDUF(γ,h,S,Ū,Z̄) + interval(1)/interval(2)*(α.^2) .* D3UF(γ,h,S,Ū,Z̄) + a[1]*D2F(m̄,γ,g,h,S,Ū,Z̄,Z̄)

    DG0 = DF(m̄,γ,g,h,S,Ū) + α .* D2F_op(m̄,γ,g,h,S,Ū,Z̄)

    component(DG0,2,2)[:,1] .= coefficients(DαG)[2:K+2]
    
    return DG0
end



function Newton_method(U, ϵ, m̄, γ, g, h, max_iter, tol)

    F_U = G(U, ϵ, m̄ , γ, g, h)
    DF_U = DG(U, ϵ, m̄ , γ, g, h)

    nf = norm(DF_U\F_U,1)
    println("Iteration 0 : ||F(U)|| = $nf")
    if nf < tol
            println("Converged in 0 iterations.")
            return U
    end
    for iter = 1:max_iter
        
        U = U - DF_U \ F_U

        F_U = G(U, ϵ, m̄ , γ, g, h)
        DF_U = DG(U, ϵ, m̄ , γ, g, h)

        nf = norm(DF_U\F_U,1)
        println("Iteration $iter: ||F(U)|| = $nf")

        if nf < tol
            println("Converged in $iter iterations.")
            return U
        end
    end
    println("Did not converge within $max_iter iterations.")
    return U
end
 




function Newton_method_bar(U, m̄, γ, g, h, max_iter, tol)

    S = space(U)
    F_U = mid.(Ḡ(U, S, m̄ , γ, g, h))
    DF_U = mid.(DḠ(U, S, m̄ , γ, g, h))

    nf = norm(DF_U\F_U,1)
    println("Iteration 0 : ||F(U)|| = $nf")
    if nf < tol
            println("Converged in 0 iterations.")
            return U
    end
    for iter = 1:max_iter
        
        U = U - DF_U \ F_U

        F_U = mid.(Ḡ(U, S, m̄ , γ, g, h))
        DF_U = mid.(DḠ(U, S, m̄ , γ, g, h))

        nf = norm(DF_U\F_U,1)
        println("Iteration $iter: ||F(U)|| = $nf")

        if nf < tol
            println("Converged in $iter iterations.")
            return U
        end
    end
    println("Did not converge within $max_iter iterations.")
    return U
end





################# Functions for the proof 






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
      







function proof_branch(K0,N,ϵ0,h,γ,m̄,g,W_fft,W_fft_big)
    
  N_ = 8N 
N_fft_ = nextpow(2, 2N_ + 1)
npts_ = N_fft_ ÷ 2 + 1

 K1 = 2*K0 + 2
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

ϵ = [interval(0.5) * ϵ0 + interval(0.5)*cospi(interval(2*k) / interval(N_fft_))*ϵ0 for k ∈ 0:npts_-1] 
ϵ = [ϵ ; reverse(ϵ)[2:npts_-1]]
m_fft_ = m̄ .- ϵ.^2


Q̄ = interval(2)*g*h .+ (m_fft_/h .- γ*h*interval(0.5)).^2
#### flat solution at the bifurcation point
ā = Sequence(Sa, interval.(zeros(K0+1))) ; ā[0] = h
Ū = [Sequence(S,[Q̄[n] ; ā[0:K0]]) for n = 1:2*(npts_-1)]

T = CartesianProduct{Tuple{ParameterSpace, CosFourier{Interval{Float64}}}}
U_fft = Vector{Sequence{T,Vector{Interval{Float64}}}}(undef, 2*(npts_-1))
J = [1:npts_-1 ; npts_+1:2*(npts_ -1)]
for n ∈ J
    U_fft[n] = Ū[n] + ϵ[n]^2*op_P(ϵ[n],S)*W_fft[n]
end
U_fft[npts_] = W_fft[npts_]


ϵb = [interval(big(0.5)) * ϵ0b + interval(big(0.5))*cospi(interval(big(2*k)) / interval(big(N_fft_)))*ϵ0b for k ∈ 0:npts_-1] 
ϵb = [ϵb ; reverse(ϵb)[2:npts_-1]]
m_fft_b = m̄b .- ϵb.^2
Q̄b = interval(big(2))*gb*hb .+ (m_fft_b/hb .- γb*hb*interval(big(0.5))).^2
#### flat solution at the bifurcation point
āb = Sequence(Sa, interval.(big.(zeros(K0+1)))) ; āb[0] = hb
Ū = [Sequence(S,[Q̄b[n] ; āb[0:K0]]) for n = 1:2*(npts_-1)]

T = CartesianProduct{Tuple{ParameterSpace, CosFourier{Interval{Float64}}}}
U_fft_big = Vector{Sequence{T,Vector{Interval{BigFloat}}}}(undef, 2*(npts_-1))
J = [1:npts_-1 ; npts_+1:2*(npts_ -1)]
for n ∈ J
    U_fft_big[n] = Ū[n] + ϵb[n]^2*op_P(ϵb[n],S)*W_fft_big[n]
end
U_fft_big[npts_] = W_fft_big[npts_]



a_fft = [Sequence(CosFourier(K0,interval(1)), component(U_fft[n],2)[:]) for n = 1:2*(npts_-1)]
Q_fft = [component(U_fft[n],1)[1] for n = 1:2*(npts_-1)]

### computation of DG on the grid

DG_grid0 = DG_grid(ϵ,h,γ,m̄,g,U_fft,CosFourier(K3, interval(1)),2*npts_-2)
#### NOTE (Sec. 7.2 / Lemma 7.3): the finite block A_K and the bound Z10 are now built
#### AFTER the tail operator A∞ below, because the new A_K depends on A∞.

#### we rename the variable n_pts_ to avoid using 2*(npts_-1) all the time
n0 = npts_
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
for n ∈ 1:npts_
    v1_grid[n] = interval(2)*(m_fft_[n]*Dy_fct_1D(e0,h) + interval(0.5)*γ*Dy_fct_1D((a_fft[n]^2),h) - γ*a_fft[n]*(Dy_fct_1D(a_fft[n],h)))
    v2_grid[n] = interval(2)*g*((Derivative(1)*a_fft[n])^2 + (Dy_fct_1D(a_fft[n],h))^2)
    v3_grid[n] =  -interval(2)*(Q_fft[n] - interval(2)*g*(a_fft[n]))*(Derivative(1)*a_fft[n])
    v4_grid[n] = -interval(2)*(Q_fft[n] - interval(2)*g*(a_fft[n]))*(Dy_fct_1D(a_fft[n],h)) 
    w[n] = mid.(im*sin2exp(v3_grid[n]) + cos2exp(v4_grid[n])) ; w[n] = real.(w[n]) ; w[n] = interval.(mid.(inv(w[n])))
    wc[n] = Sequence(CosFourier(K0,interval(1)),real.((interval.(0.5)*(w[n] + reverse(w[n])))[0:K0]))
    ws[n] = Sequence(SinFourier(K0,interval(1)),real.((interval.(0.5)*(w[n] - reverse(w[n])))[1:K0]))
end

#### inverse of the derivative (will be used for high frequencies)
Dx_inv = LinearOperator(SinFourier(K4,interval(1)), Sa4, interval.(zeros(K4+1,K4)))
for k = 1:K4
    Dx_inv[k,k] = - interval(1)/interval(k)
end
### operator A∞ = π^{>K1} D_y^{-1} M_w on the grid (tail rows > K1).
### NEW: now computed for ALL n (including the bifurcation index n0), since the tail A∞
### also enters the construction of the finite block A_K below (Sec. 7.2, eq. 7.14).
Dy_Mw_grid =  Vector{LinearOperator{CosFourier{Interval{Float64}},CosFourier{Interval{Float64}},Matrix{Interval{Float64}}}}(undef, npts_)
for n ∈ 1:npts_
    Dy_Mw_grid[n] = Dy_op_inv(Sa4,h)*project(Multiplication(wc[n]),Sa3,Sa4) - Dx_inv*project(Multiplication(ws[n]),Sa3,SinFourier(K4,interval(1)))
    Dy_Mw_grid[n][0:K1,:] .= interval.(zeros(K1+1,K3+1)) ### we only keep the tail of the operator
end

###################################################################################
#### NEW (Sec. 7.2, eq. 7.14): the finite block A_K is the numerical inverse of  ##
####     M = π^{≤K1} DG (I − A∞ DG) π^{≤K1}.                                      ##
#### Only the a-row of M receives the A∞ correction: the Q-row is unchanged       ##
#### because A∞ DG is supported on modes > K1 while the Q-row only reads mode 0.   ##
###################################################################################
M_grid0 = [ begin
        DGaa = component(DG_grid0[n],2,2)                                     # Sa2 -> Sa2
        DGa1 = component(DG_grid0[n],2,1)                                     # ParameterSpace -> Sa2
        AinfDGaa = Dy_Mw_grid[n]*project(DGaa,Sa1,Sa2)                        # A∞ DG (a-col) : Sa1 -> Sa4
        AinfDGa1 = Dy_Mw_grid[n]*project(DGa1,ParameterSpace(),Sa2)           # A∞ DG (Q-col) : Param -> Sa4
        Cor_aa = project(DGaa,Sa2,Sa1)*project(AinfDGaa,Sa1,Sa2)             # π DG A∞ DG (a->a) : Sa1 -> Sa1
        Cor_a1 = project(DGaa,Sa2,Sa1)*project(AinfDGa1,ParameterSpace(),Sa2) # π DG A∞ DG (Q->a) : Param -> Sa1
        M = project(DG_grid0[n],S1,S1)                                        # π^{≤K1} DG π^{≤K1}
        component(M,2,2) .= coefficients(component(M,2,2)) - coefficients(Cor_aa)
        component(M,2,1) .= coefficients(component(M,2,1)) - coefficients(Cor_a1)
        M
    end for n = 1:npts_ ]

#### A_K = numerical inverse of M, re-expanded as a Chebyshev polynomial of order N (DFT)
A_grid = inv.([mid.(M_grid0[i]) for i=1:npts_])

A_cheb = [interval.(mid.(rifft!(complex.([A_grid[k][i,j] for k = 1:length(A_grid)]), Chebyshev(N)))) for i = indices(codomain(A_grid[1])), j = indices(domain(A_grid[1]))]

A_grid = real.(cheb2grid(A_cheb, N_fft_))
A_grid = [LinearOperator(S1,S1 ,A_grid[i]) for i = 1:npts_]

######### Computation of the Z1 bound ###########

ADG_grid = A_grid.*DG_grid0

#### NEW Z10 (Lemma 7.3): Z10 = || π^{≤K1} − A_K P DF (π^{≤K1} − A∞ DF) P π^{≤K1} ||_H1
####                           = || I − A_K M ||,   M = π^{≤K1} DG (I − A∞ DG) π^{≤K1}
AM_grid = A_grid.*M_grid0
Id = LinearOperator(S1, S1, interval.(1.0*I[1:K1+2, 1:K1+2]))
ADG_Z_grid = [Id-D0_op(S1)*project(AM_grid[i],S1,S1)*D0_inv_op(S1) for i = 1:npts_]

Z10 = norm_op_grid(ADG_Z_grid, N_, 1, 1)

println("Bound Z10 : Z10 = $Z10")

#### we use the tail operator A∞ for computing the bound Z11. We initialize Dy^{-1}*Mw*DG
CDG_grid = Vector{LinearOperator{CosFourier{Interval{Float64}},CosFourier{Interval{Float64}},Matrix{Interval{Float64}}}}(undef, npts_)
### computation of Dy^{-1}*Mw*DG on the grid. We also incorporate the weights in the norm given by D0 and D0_inv
for n ∈ J
    CDG_grid[n] = D0_op(Sa4)*Dy_Mw_grid[n]*project(component(DG_grid0[n],2,2),Sa1,Sa3)*D0_inv_op(Sa1)
end
CDG_grid[n0] = interval(0)*CDG_grid[n0+1]
Z11 = norm_op_grid(CDG_grid, N_, K1+2, 1)

println("Bound Z11 : Z11 = $Z11")

#### NEW Z12 (Lemma 7.3): Z12 = || A_K(ε) DF(V,m̄−ε²) π^{>K1} ||_H1   (this is the former Z14)
ADG_Z_grid = [D0_op(S1)*ADG_grid[i]*D0_inv_op(S3) for i = 1:npts_]
ADG_Z_grid[n0] = interval(0)*ADG_Z_grid[n0+1]
Z12 = norm_op_grid(ADG_Z_grid, N_, 1, K1+2)
println("Bound Z12 : Z12 = $Z12")


##### Computation of the different components of Z∞
### We start by computing δ
w_prod = [-ws[n]*v3_grid[n] + wc[n]*v4_grid[n] - interval(1) for n ∈ 1:npts_]
w_v4 = [w[n]*cos2exp(v4_grid[n]) for n ∈ 1:npts_]

#### NEW (Lemma 7.3): the exponential decay rate is now e^{-2 d K1} (was e^{-2 d K0})
δ = norm_seq_grid(w_prod, N_) +  interval(2)*norm_seq_grid(w_v4, N_)*exp(-interval(2*K1)*h)/(interval(1) - exp(-interval(2*K1)*h))
δ = δ*interval(K1+2)/interval(K1+1) 


D = LinearOperator(Fourier(K0,interval(1)),Fourier(K0,interval(1)), interval.(1.0*I[1:2K0+1, 1:2K0+1]))
for k = -K0:-1
    D[k,k] = interval(k) - abs(interval(k))*coth(h*abs(interval(k)))
end
for k = 1:K0
    D[k,k] = interval(k) - abs(interval(k))*coth(h*abs(interval(k)))
end
    D[0,0] = - interval(1)/h

#### first term of Z∞   (NEW coefficient 1/(K1+1), was (K1+2)/(K1+1)^2 — Lemma 7.3)
w_prod = [w[n]*(cos2exp(v1_grid[n])*(D*cos2exp(a_fft[n])) + cos2exp(v2_grid[n])) for n ∈ 1:npts_]
w_prod[n0] = interval(0)*w_prod[n0+1]
Z∞ = norm_seq_grid(w_prod,N_)/interval(K1+1)
### adding the second term
w_prod = [w[n]*cos2exp(v1_grid[n]) for n ∈ 1:npts_]
w_prod[n0] = interval(0)*w_prod[n0+1]
a_fft_D = [D0_op(Sa1)*a_fft[n] for n ∈ 1:npts_]
a_fft_D[n0] = interval(0)*a_fft_D[n0+1]

Z∞ = Z∞ + interval(4)*norm_seq_grid(w_prod,N_)*exp(-interval(2*K0)*h)/(interval(1) - exp(-interval(2*K0)*h))*interval(K0)*interval(K1+2)/interval(K1+1)*norm_seq_grid(a_fft_D,N_)
Z∞ = Z∞ + δ #### we finish by adding epsilon computed above

println("Bound Z∞ : Z∞ = $Z∞")

#### NEW Z1 (eq. 7.16): Z1 = max{ Z10 + Z11 , (1 + Z12) Z∞ }
Z1 = maximum([Z10+Z11 (interval(1)+Z12)*Z∞])

println("Bound Z1 : Z1 = $Z1")


########## Computation of the Y bound ###########       
#### computation of AK*F(U,m) on the grid of m's
G_grid0, F_raw =  G_grid(ϵb,hb,γb,m̄b,gb,U_fft_big,S4,npts_)
#### Y FIRST TERM (Lemma 7.4): (1/eps^2)||A_K (I - DG A∞) P F||_H1 = ||A_K (I - DG A∞) G||_H1.
#### Only the a-row receives the (DG A∞) correction (A∞ G is a tail, so the Q-row reads mode 0 -> 0),
#### exactly as in the construction of A_K / M used for Z10.

### recompute the action of A∞ for higher frequencies argument (since the input is of order K4)
Dy_Mw_grid =  Vector{LinearOperator{CosFourier{Interval{Float64}},CosFourier{Interval{Float64}},Matrix{Interval{Float64}}}}(undef, npts_)
for n ∈ 1:npts_
    Dy_Mw_grid[n] = Dy_op_inv(Sa4,h)*project(Multiplication(wc[n]),Sa4,Sa4) - Dx_inv*project(Multiplication(ws[n]),Sa4,SinFourier(K4,interval(1)))
    Dy_Mw_grid[n][0:K1,:] .= interval.(zeros(K1+1,K4+1)) ### we only keep the tail of the operator
end

f_a = [component(F_raw[i],2) for i = 1:npts_]                           # raw residual f = F(V,mbar-eps^2)
AGG = [ begin
        aG      = component(G_grid0[n],2)                              # Sa4   (a-component of G)
        AinfG   = Dy_Mw_grid[n]*aG                        # A∞ aG  (tail, Sa4)
        DGAinfG = project(component(DG_grid0[n],2,2),Sa3,Sa1)*project(AinfG,Sa3)   # DG[2,2] A∞ aG : Sa1
        GG = project(G_grid0[n],S1)                                    # pi^{<=K1} G
        component(GG,2) .= coefficients(component(GG,2)) - coefficients(DGAinfG)
        A_grid[n]*GG                                                   # A_K (I - DG A∞) G  on S1
    end for n = 1:npts_ ]
Y = norm_seq_grid([D0_fct(AGG[i]) for i = 1:npts_],N_)

#### Y SECOND TERM (Lemma 7.4): (K1+2)/(K1+1) ( ||pi^{k>K1}(w*f)|| + ||pi^{k<-K1}(conj(w)*f)|| ).
#### w has real coefficients, so the conjugate w* corresponds to reverse(w).
wf_pos = [w[i]*cos2exp(f_a[i]) for i = 1:npts_]                         # w * f
wf_neg = [reverse(w[i])*cos2exp(f_a[i]) for i = 1:npts_]                # conj(w) * f
#### keep only the requested one-sided tail. Fourier coefficients are ordered mode -M .. M,
#### so mode k sits at vector index k+M+1.
for i = 1:npts_
    Mp = order(wf_pos[i]) ; coefficients(wf_pos[i])[1:(K1+Mp+1)]        .= interval(0)   # zero k <= K1  (keep k > K1)
    Mn = order(wf_neg[i]) ; coefficients(wf_neg[i])[(Mn-K1+1):(2*Mn+1)] .= interval(0)   # zero k >= -K1 (keep k < -K1)
end
Y = Y + interval(K1+2)/interval(K1+1)*( norm_seq_grid(wf_pos,N_) + norm_seq_grid(wf_neg,N_) )

println("Bound Y : Y = $Y")

######## computation Z2 bound (Lemma 7.5) ##############

r0 = interval(1e-4)

norm_ell_1_a = norm_seq_grid(a_fft,N_)
norm_H1_a = norm_seq_grid(a_fft_D,N_)

#### NEW Z_D2_1 (Lemma 7.5): Z_D2_1 = sup || A P D2F(V,mbar-eps^2)(Zbar,Zbar) ||_H1.
#### On span(Zbar): P D2F(Zbar,Zbar) = D2F(Zbar,Zbar) and A acts as A_K (the tail Ainf
#### vanishes on these low modes), so we apply the current A_grid (= A_K, before rescaling).
Zb4 = Sequence(S4, interval.(zeros(K4+2))) ; component(Zb4,2)[1] .= interval(1)
AD2 = [A_grid[n]*project(D2F(m_fft_[n],γ,g,h,S4,project(U_fft[n],S4),Zb4,Zb4),S1) for n = 1:npts_]
ZD21 = norm_seq_grid([D0_fct(AD2[i]) for i = 1:npts_], N_)
println("Bound ZD2_1 : ZD2_1 = $ZD21")

####### computation of Z_A = sup eps||W A P||  (norm of the approximate inverse) #########

for n = J
    A_grid[n] = D0_op(S1)*A_grid[n]*op_P(ϵ[n],S1)*ϵ[n]
end
component(A_grid[n0],2,2)[:,[0 ; 2:K1]] = interval.(zeros(K1+1,K1)) 
component(A_grid[n0],1,2)[1,[0 ; 2:K1]] = interval.(zeros(1,K1)) 

norm_A = norm_op_grid(A_grid, N_, 1, 1)
norm_w = norm_seq_grid(w,N_)

norm_A = norm_A + (I"1" + Z12)*interval(K1+1)/interval(K1+1)*norm_w ## overestimation of the norm of A, but it is enough for the proof

######## computation of ZD2 = Z_D2_2 ##############

fct1 = [m_fft_[n]/h + interval(0.5)*γ*Dy_fct_1D((a_fft[n]^2),h) - γ*a_fft[n]*(Dy_fct_1D(a_fft[n],h)) for n = 1:npts_]
fct2 = [Q_fft[n]-interval(2)*g*(a_fft[n]) for n = 1:npts_]

ZD2 = interval(2)*γ^2*(interval(2)*norm_H1_a + norm_ell_1_a)^2/h^2 
ZD2 += interval(6)*abs(γ)/h*norm_seq_grid(fct1,N_)
#### NEW (Lemma 7.5): factor 2 in front of ||Qbar e0 - 2g abar||
ZD2 += (interval(1) + interval(1)/h^2)*(interval(4)*(interval(1)+interval(2)*g)*norm_H1_a + interval(2)*norm_seq_grid(fct2,N_))

###### computation of ZD3 = Z_D3(r0)  (the leading r0 factor is now applied in Z2 below) #############
max_Q = norm(rifft!(complex.(Q_fft), Chebyshev(N_)),1) 

ZD3 = interval(27)*γ^2/h^2*(norm_H1_a+r0) + interval(4)*g*(interval(1) + interval(1)/h^2) 
ZD3 += (abs(max_Q) +interval(2)*g)*(interval(1) + interval(1)/h^2) 
ZD3 += interval(3)*(abs(max_Q)+r0)*(interval(1) + interval(1)/h^2)

#### NEW Z2 (Lemma 7.5): Z2(r0) = max{ Z_D2_1, Z_A Z_D2_2, eps0 Z_A Z_D2_2 } + Z_A Z_D3(r0) r0
Z2 = maximum([ZD21  norm_A*ZD2  ϵ0*norm_A*ZD2]) + norm_A*ZD3*r0
println("Bound Z2 : Z2 = $Z2")


if inf((interval(1)-Z1)^2) > sup(interval(2)*Z2*Y)
    rm = (interval(1)-interval(sup(Z1)) - sqrt((interval(1)-interval(sup(Z1)))^2-interval(2)*Z2*Y))/Z2
    rp = (interval(1)-interval(sup(Z1)) + sqrt((interval(1)-interval(sup(Z1)))^2-interval(2)*Z2*Y))/Z2
    rp = minimum([rp (interval(1)-Z1)/Z2 r0])
    if (sup(abs(rm))<inf(r0))&&(sup(Z1 + Z2*rm) < 1)
        println("Proof successful with radius in : [$rm , $rp ]")
        return rm,rp
    else
        display("second condition not verified")
        return naN
    end
else
    display("first condition not verified")
    return naN
end  

end


cheb_fourier(grid, Sx, N, npts_) =
    Sequence(Chebyshev(N) ⊗ Sx, reduce(vcat, [coefficients(real.(rifft!(complex.([coefficients(grid[j])[k] for j = 1:2*(npts_-1)]), Chebyshev(N)))) for k = 1:length(coefficients(grid[1]))]))


function proof_geometry(K0,N,ϵ0,h,γ,m̄,g,W_fft,rm)


    ϵ = [interval(0.5) * ϵ0 + interval(0.5)*cospi(interval(2*k) / interval(N_fft_))*ϵ0 for k ∈ 0:npts_-1] 
    ϵ = [ϵ ; reverse(ϵ)[2:npts_-1]]
    m_fft_ = m̄ .- ϵ.^2


    Q̄ = interval(2)*g*h .+ (m_fft_/h .- γ*h*interval(0.5)).^2
    #### flat solution at the bifurcation point
    ā = Sequence(Sa, interval.(zeros(K0+1))) ; ā[0] = h
    S = ParameterSpace()×CosFourier(K0,interval(1)) 
    Ū = [Sequence(S,[Q̄[n] ; ā[0:K0]]) for n = 1:2*(npts_-1)]

    T = CartesianProduct{Tuple{ParameterSpace, CosFourier{Interval{Float64}}}}
    U_fft = Vector{Sequence{T,Vector{Interval{Float64}}}}(undef, 2*(npts_-1))
    J = [1:npts_-1 ; npts_+1:2*(npts_ -1)]
    for n ∈ J
        U_fft[n] = Ū[n] + ϵ[n]^2*op_P(ϵ[n],S)*W_fft[n]
    end
    U_fft[npts_] = W_fft[npts_]

    a_fft = [Sequence(CosFourier(K0,interval(1)), component(U_fft[n],2)[:]) for n = 1:2*(npts_-1)]

    ξx = [Dy_fct_1D(a_fft[n],h) for n = 1:2*(npts_-1)]
    ξx = cheb_fourier(ξx,CosFourier(K0,interval(1)),N,npts_)

    dX = 0.01
    dE = 0.01
    nx = ceil(Int, π/dX) ; xnodes = range(-π, 0, length = nx+1)        # consecutive nodes tile (-π,0)
    ne = ceil(Int, 2/dE)  ; enodes = range(-1, 1, length = ne+1)       # consecutive nodes tile [-1,1] (the Chebyshev variable ↔ ϵ ∈ [0,ϵ0])
    S = interval(1)
    for j = 1:ne
        eI = interval(enodes[j], enodes[j+1])
        xI = interval(inf(interval(-π)),sup(interval(-π)) + sup(I"1.1"*dX))   # handle the case -π separately
        S = minimum([S ξx(eI,xI)-coth(h)*rm])               
        for i = 2:nx
            xI = interval(xnodes[i], xnodes[i+1]) 
            S = minimum([S ξx(eI,xI)-coth(h)*rm])
        end
    end

    if inf(S) > 0
        println("The branch satisfies ξx >0 on (-π,0) for all ϵ ∈ [0,ϵ0]")
    else 
        println("Cannot conclude. Refine grid. value of S = $S")
    end 


end 
