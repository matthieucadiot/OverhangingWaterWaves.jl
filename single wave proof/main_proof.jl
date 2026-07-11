using RadiiPolynomial, IntervalArithmetic, LinearAlgebra, Plots, JLD2
using Logging ; disable_logging(Logging.Info)

include("list_fct_proof.jl")

K0 = 200 ; K1 = 2*K0 + 2 ; γ = -interval(5) ; g = interval(1) ; h = interval(2) 
setprecision(256)

U = load("U.jld2","U_big")
m_big = - parse(Interval{BigFloat}, "0.85")
m = -I"0.85"

U_big = U  #### the approximate solution is saved in bigfloats
U = interval.(Float64.(inf.(U),RoundDown),Float64.(sup.(U),RoundUp) ) ### we save a copy in float for computations where the bigfloats are not needed

### various spaces of sequences with different sizes for the Fourier modes
Sa = CosFourier(K0,interval(1))
Sa1 = CosFourier(K1,interval(1))
Sa2 = CosFourier(K1+2*K0,interval(1))
Sa3 = CosFourier(K1+3*K0,interval(1))
Sa4 = CosFourier(K1+4*K0,interval(1))
### corresponding product spaces with a parameter Q
S = ParameterSpace()×Sa
S1 = ParameterSpace()×Sa1
S2 = ParameterSpace()×Sa2
S3 = ParameterSpace()×Sa3
S4 = ParameterSpace()×Sa4


### components of U
Q = component(U,1)[1]
a = component(U,2) ; a = Sequence(Sa,vec(coefficients(a)))
### components of U_big, only a_big is needed. Q does not have to be in big floats
a_big = component(U_big,2)
Q_big = component(U_big,1)[1]
#### construction of the linear operators used in the proof, all of them are diagonal in the Fourier basis
K2 = K1 + 2*K0 ; K3 = K1 + 3*K0 ; K4 = K1 + 4*K0;
n = interval.(1:K3) ; n4 = interval.(1:K4)
Id = LinearOperator(S1, S1, interval.(1.0*I[1:K1+2, 1:K1+2]))
Dy = LinearOperator(Sa3, Sa3, Diagonal( vec([interval(1)/h ; n.*coth.(h*n)]) ))
Dy_inv_hat = LinearOperator(Sa4, Sa4, Diagonal( vec([h ; interval(1) ./ n4]) ))
#### now we define weights for the norms on X
Ŵ = LinearOperator(Sa3, Sa3, Diagonal( vec([interval(1) ; interval(2)*n]) ))
Ŵ_inv = LinearOperator(Sa3, Sa3, Diagonal( vec([interval(1) ; interval(0.5) ./ n]) ))
### weights for the norms on H_1
𝒲 = LinearOperator(S3, S3, Diagonal( vec([interval(1) ; interval(1) ; interval(2)*n]) ))
𝒲_inv = LinearOperator(S3, S3, Diagonal( vec([interval(1) ; interval(1) ; interval(0.5) ./ n]) ))

Dx_inv = LinearOperator(SinFourier(K4,interval(1)), Sa4, interval.(zeros(K4+1,K4)))
for k = 1:K4
    Dx_inv[k,k] = - interval(1)/interval(k)
end
Dx = project(Derivative(1), Sa3, SinFourier(K3,interval(1)))



##### computation of approximate inverse ########

### construction of the functions v_i
e0 = Sequence(Sa, interval.(zeros(size(a)[1])))
e0[0] = interval(1)
v1 = interval(2)*(m*Dy*e0 + interval(0.5)*γ*Dy*(a^2) - γ*a*(Dy*a))
v2 = interval(2)*g*((Dx*a)^2 + (Dy*a)^2)
v3 = -interval(2)*(Q - interval(2)*g*(a))*(Dx*a)
v4 = -interval(2)*(Q - interval(2)*g*(a))*(Dy*a) 

v1 = project(v1,CosFourier(2*K0,interval(1))) ; v2 = project(v2,CosFourier(2*K0,interval(1))) ;v3 = project(v3,SinFourier(2*K0,interval(1))) ;v4 = project(v4,CosFourier(2*K0,interval(1)))

### construction of the function w, approximating the inverse of the function v4 + iv3
w = mid.(im*sin2exp(v3) + cos2exp(v4)) ; w = real.(w) ; w = inv(w) ; 
w = interval.(w) 
### we recover the cosine and the sine parts of w. Since the proof uses cosine series only, this decommposition will be useful in computations
wc = interval(0.5)*(w + reverse(w)) ; ws = interval(0.5)*(w - reverse(w))
wc = Sequence(CosFourier(2*K0,interval(1)),wc[0:end]) ; ws = Sequence(SinFourier(2*K0,interval(1)), ws[1:2*K0])

### operator C
Ma = project(Multiplication(a), Sa3, Sa3)
C = γ*(Dy*Ma - Ma*Dy - project(Multiplication(Dy*a), Sa3, Sa3))
### jacobian Daf
Daf0 = project(Multiplication(v1), Sa3, Sa3)*C + project(Multiplication(v2), Sa3, Sa3) + project(Multiplication(v3), SinFourier(K1+3*K0,interval(1)), Sa3)*Dx + project(Multiplication(v4), Sa3, Sa3)*Dy
DQf0 = - project(v2,Sa3)/(interval(2)*g) ### partial derivative with Q

#### computation of the jacobian of F
DF0 = LinearOperator(S3, S3, interval.(zeros(K3+2, K3+2)))
component(DF0,2,2).= coefficients(Daf0)
component(DF0,2,1).= coefficients(project(DQf0,Sa3))
component(DF0,1,2).= coefficients(project(e0,Sa3))'

##### construction of the approximate inverse  A = A_K - A_K DF(U) A∞ + A∞  (Schur complement, eq (6.34)) #####

### diagonal operator 𝒟,  (𝒟a)_k = k - ω_k  (Lemma 6.3), used in the Z∞ bound
n0 = interval.(1:K0)
𝒟 = LinearOperator(Sa, Sa, Diagonal( vec([-interval(1)/h ; n0 - n0.*coth.(h*n0)]) ))

### analytic tail inverse  Â∞ = π^{>K1} D̂y⁻¹ M_w  (eq (6.29)), restricted to the tail K1 < |k| ≤ K1+2K0
Â∞ = Dy_inv_hat*project(Multiplication(wc),Sa3,Sa3) - Dx_inv*project(Multiplication(ws),Sa3,SinFourier(K1+3*K0,interval(1)))
Â∞ = project(Â∞,Sa3,Sa3) ; Â∞[0:K1,:] .= interval.(zeros(K1+1,K3+1))
### embedding of the tail inverse into H1 :  A∞ U = (0, Â∞ a)  (eq (6.33))
A∞ = LinearOperator(S3, S3, interval.(zeros(K3+2, K3+2)))
component(A∞,2,2) .= coefficients(Â∞)

### Schur-complement operator  M = DF(U) - DF(U) A∞ DF(U)
M = DF0 - DF0*A∞*DF0

##### finite dimensional truncation for the inverse :  A_K = π^{≤K1} A_K π^{≤K1}
A_K = interval.(inv(mid.(project(M,S1,S1))))

### the components of A_K, A∞ are now computed and we are ready to compute the bounds for the radii polynomial
########### Computation of the bounds for the radii polynomial ############

println("\n" * "="^74)
println("CONSTRUCTIVE EXISTENCE PROOF of a single overhanging wave  (m = $(mid(m)))")
println("  Computing the Newton–Kantorovich bounds (Y, Z₁, Z₂) ...")
println("-"^74)


##### computation of Z1 

AKM = A_K*M
Z10 = opnorm( Id - 𝒲*project(AKM,S1,S1)*𝒲_inv, 1)
println("Bound Z10 : Z10 = $Z10")

𝒞 = project(Ŵ*Â∞*project(Daf0,Sa1,Sa3)*Ŵ_inv,Sa1,Sa4)
Z11 = opnorm(𝒞[K1+1:end,:],1)
println("Bound Z11 : Z11 = $Z11")

##### Z12 = ‖A_K DF(U) π^{>K1}‖_{H1}  (eq (6.35), plays the role of the former Z14) #####
AKDF = A_K*DF0
Z12 = opnorm(component(𝒲*AKDF,2,2)[0:K1,K1+1:K1+2*K0].*vec(interval(0.5) ./interval(K1+1:K1+2*K0))',1)
println("Bound Z12 : Z12 = $Z12")

δ = interval(K1+2)/interval(K1+1)*(norm(-ws*v3 + wc*v4 - interval(1),1) + norm(w*cos2exp(v4),1)*interval(2)*exp(-interval(2)*interval(K1)*h)/(interval(1)-exp(-interval(2)*interval(K1)*h)))

Z∞ = δ + norm(w*(cos2exp(v1*(𝒟*a)) + cos2exp(v2)), 1)/interval(K1+1) + interval(K1+2)/interval(K1+1)*norm(w*cos2exp(v1),1)*interval(4)*interval(K0)*exp(-interval(2)*interval(K0)*h)/(interval(1)-exp(-interval(2)*interval(K0)*h))*norm(Ŵ*a,1) 

println("Bound Z∞ : Z∞ = $Z∞")

Z1 = maximum([Z10+Z11 (interval(1)+Z12)*Z∞])
println("Bound Z1 : Z1 = $Z1")

###### computation of the Y bound 
n = interval.(big.(1:K2))
Dy_big = LinearOperator(Sa2, Sa2, Diagonal( vec([I"1"/h ; n.*coth.(h*n)]) ))

f0 = (m_big*project(Dy*e0,Sa) + I"0.5"*γ*project(Dy_big*(a_big^2),CosFourier(2*K0,interval(1))) - γ*a_big*(project(Dy_big*a_big,Sa)))^2 - (Q_big[1] - I"2"*g*a_big)*((project(Derivative(1)*a_big,SinFourier(K0,interval(1)))^2 + (project(Dy_big*a_big,Sa))^2))

F0 = Sequence(ParameterSpace()×space(f0), [a_big[0]-h; coefficients(f0)])

F0_S2 = project(F0,S4)
Y1 = A_K*(F0_S2 - DF0*(A∞*F0_S2))                  # A_K (I - DF(U) A∞) F(U)
wf = w*cos2exp(f0) ; wfc = reverse(w)*cos2exp(f0) ; Nw = order(wf)   # w*f(U)  and  w⋆*f(U)
Y = norm(𝒲*Y1,1) + interval(K1+2)/interval(K1+1)*( norm(wf[K1+1:Nw],1) + norm(wfc[-Nw:-(K1+1)],1) )
println("Bound Y : Y = $Y")


######## computation of the Z2 bound

### we suppose a priori that the radius of contraction will be smaller than r_Z2 = 1e-6. This will have to be verified a posteriori, but it allows us to simplify some of the bounds by neglecting higher order terms in r.
r_Z2 = interval(1e-6)

max_h = maximum([interval(1) interval(1)/h])

A = project(A_K,S1,S1) - project(A_K*DF0*A∞,S1,S1) + project(A∞,S1,S2)
norm_A = maximum([opnorm(𝒲*A,1)  (interval(1)+Z12)*interval(K1+2)/interval(K1+1)*norm(w,1)])

ZD2 = interval(2)*γ^2*max_h^2*(interval(2)*norm(Ŵ*a,1) + norm(a,1))^2 
ZD2 += interval(6)*max_h*abs(γ)*norm(m*Dy*e0 + interval(0.5)*γ*Dy*(a^2) - γ*a*(Dy*a),1) 
ZD2 += (interval(1) + max_h^2)*(interval(4)*(interval(1)+interval(2)*g)*norm(Ŵ*a,1) + interval(2)*norm(Q*e0 - I"2"*g*a,1) )


ZD3 = interval(27)*γ^2/h^2*(norm(Ŵ*a,1) + r_Z2) + interval(4)*g*(interval(1) + max_h^2)
ZD3 += (abs(Q)+interval(2)*g)*(interval(1) + max_h^2) 
ZD3 += interval(3)*(abs(Q)+r_Z2)*(interval(1) + max_h^2)

Z2 = norm_A*(ZD2 + ZD3*r_Z2)
println("Bound Z2 : Z2 = $Z2")

##### Radii polynomial and conclusion #####

r0 = (1 - sup(big(Z1)) - sqrt( (1-sup(big(Z1)))^2 - 2*sup(big(Z2))*sup(Y) ) )/inf(Z2)
r0 = interval(r0)

p1 = I"0.5"*Z2*r0^2 - (I"1" - Z1)*r0 + Y
p2 = Z1 + Z2*r0

if sup(p1) < 0

    if (sup(p2) < 1)&&(sup(r0)<inf(r_Z2))
        sup_r = Float64(sup(r0),RoundUp)
        println("-"^74)
        println("  ✓ EXISTENCE PROVEN : a unique wave exists near the approximate solution.")
        println("    Newton–Kantorovich bounds :  Y  = $(sup(Y))")
        println("                                 Z₁ = $(sup(Z1))")
        println("                                 Z₂ = $(sup(Z2))")
        println("    Validated radius          :  r₀ = $sup_r   (contraction radius / proven accuracy)")
        println("    ⇒ the exact wave lies within distance r₀ of the approximation in the norm of X.")
    else
        println("  ✗ second condition (p₂ < 1 and r₀ < r_Z2) not verified — proof inconclusive.")
        return naN
    end
else
    println("  ✗ first condition (radii polynomial p₁ < 0) not verified — proof inconclusive.")
    return naN
end



println("\n" * "="^74)
println("VERIFYING THE GEOMETRY of the proven wave  (η>0, monotonicity, overhang, injectivity)")
println("-"^74)
###### We prove that the profile satisfies the non-degeneracy conditions
#### First we prove that η>0
dx = 0.01
X = -π:dx:π

S = interval(1)
for x ∈ X 
    global S
    xI = interval(x-dx,x+dx) 
    S = minimum([S a(xI)-r0])  ## add the error of the approximate solution r0
end 

if inf(S) > 0
    println("The profile satisfies η>0")
else 
    println("Could not conclude if η>0")
end 

#### We prove monotonicity on (-π,0)
### Since ηx(-π) = ηx(0), we need to control the second order derivative ηxx around -π and around 0

dx = 0.001
ϵ = 0.1
X = -π+ϵ:dx:-ϵ  

## verification of positivity on a smaller interval, given by ϵ
S = interval(1)
for x ∈ X 
    global S
    xI = interval(x-dx,x+dx) 
    S = minimum([S (Derivative(1)*a)(xI)-r0]) ## add the error of the approximate solution r0
end 


if inf(S) > 0
    println("The profile satisfies ηx>0 on X = (-π + ϵ, -ϵ) with ϵ = $ϵ")
else 
    println("Could not conclude if ηx>0 on X = (-π + ϵ, -ϵ) with ϵ = $ϵ")
end 


### we prove that ηxx is strictly positive on (-π,-π+ϵ) and on (-ϵ,0). This will imply that ηx is monotonely increaing on (-π,0). For this we use the function bound_error_axx(), which computes explicitly the error bound in C^0 between D_x^2 ã and its approximation given in Lemma 7.12

āxx =  axx(a, Q, m, h, γ, g, w)
err = bound_error_axx(a, Q, m, h, γ, g, w, δ, r0) ### The error is uniform on [-π,π] and computed thanks to Lemma 7.12

dx = 0.001
X = -π:dx:-π+ϵ 
S = interval(1)
for x ∈ X 
    global S, err
    xI = interval(x-dx,x+dx) 
    S = minimum([S āxx(xI)-err]) ## add the error of the approximate second derivative
end 


if inf(S) > 0
    println("The profile satisfies ηxx>0 on X = (-π ,-π + ϵ) with ϵ = $ϵ")
else 
    println("Could not conclude if ηxx>0 on X = (-π ,-π + ϵ) with ϵ = $ϵ")
end 



dx = 0.01
X = -ϵ:dx:0
S = interval(-1)
for x ∈ X 
    global S, err
    xI = interval(x-dx,x+dx) 
    err = bound_error_axx(a, Q, m, h, γ, g, w, δ, r0)
    S = maximum([S āxx(xI)+err]) ## add the error of the approximate second derivative
end 


if sup(S) < 0
    println("The profile satisfies ηxx < 0 on X = (-ϵ , 0) with ϵ = $ϵ")
else 
    println("Could not conclude if ηxx < 0 on X = (-ϵ , 0) with ϵ = $ϵ")
end 


##### We prove that the wave satisfies -π < ξ(x) < 0 for all x ∈ (-π, 0) in order to prove injectivity of the profile. We start by verifying that  -π < ξ(x) < 0 for all x ∈ X = (-π + ϵ, -ϵ), and then prove that ξx is positive at (-π, -π + ϵ) and positive at (-ϵ, 0). Using that ξ(-π) = -π and ξ(0) = 0, this provides the desired proof.

Hhinv = LinearOperator(Sa, SinFourier(K0,interval(1)), interval.(zeros(K0,K0+1)) )
for k=1:K0
    Hhinv[k,k] = coth(interval(k)*h) ## Hilbert transform from cosine to sine series
end

ξ_minus_x = Hhinv*a ## ξ = x + Hhinv*η, we compute only the part Hhinv*η here

dx = 0.001
ϵ = 0.1
X = -π+ϵ:dx:-ϵ  

## verification of positivity on a smaller interval, given by ϵ
S1 = interval(1)
S2 = interval(-1)
for x ∈ X 
    global S1, S2
    xI = interval(x-dx,x+dx) 
    S1 = minimum([S1 interval(π)+xI+ξ_minus_x(xI)-coth(h)*r0]) ## add the error of the approximate solution r0
    S2 = maximum([S2 xI+ξ_minus_x(xI)+coth(h)*r0])
end 


if (inf(S1) > 0)&&(sup(S2) < 0)
    println("The profile satisfies -π < ξ(x) < 0 on X = (-π + ϵ, -ϵ) with ϵ = $ϵ")
else 
    println("Could not conclude if -π < ξ(x) < 0 on X = (-π + ϵ, -ϵ) with ϵ = $ϵ")
end 



### verification that ξx = Dy*η is positive at (-π, -π + ϵ) and negative at (-ϵ, 0)

ξx = Dy_fct_1D(a,h)
dx = 0.01
X = -π:dx:-π+ϵ 
S = interval(1)
for x ∈ X 
    global S
    xI = interval(x-dx,x+dx) 
    S = minimum([S ξx(xI)-coth(h)*r0]) 
end 


if inf(S) > 0
    println("The profile satisfies ξx>0 on X = (-π ,-π + ϵ) with ϵ = $ϵ")
else 
    println("Could not conclude if ξx>0 on X = (-π ,-π + ϵ) with ϵ = $ϵ")
end 


dx = 0.01
X = -ϵ:dx:0
S = interval(1)
for x ∈ X
    global S 
    xI = interval(x-dx,x+dx) 
    S = minimum([S ξx(xI)-r0]) 
end 


if inf(S) > 0
    println("The profile satisfies ξx > 0 on X = (-ϵ , 0) with ϵ = $ϵ")
else 
    println("Could not conclude if ξx > 0 on X = (-ϵ , 0) with ϵ = $ϵ")
end 



#### Finally, we prove that the profile is overhanging. For this, we prove that there exists some x ∈ (-π,0) such that ξx(x) < 0

X = -π:dx:0  
for x ∈ X 
    xI = interval(x,x+dx) 
    if sup(ξx(xI) + coth(h)*r0) < 0
        return println("The wave is overhanging.")
    end
end 



#### we plot the wave profile on 2 periods and save as a pdf
a0 = mid.(a) ; a0 = Sequence(CosFourier(K0,1),coefficients(a0)) ; m = mid(m)

X = -2π:0.01:2π
Na = size(a0)[1] - 1
Cd = LinearOperator(CosFourier(Na,1) , SinFourier(Na,1), [zeros(Na) Matrix(Diagonal( coth.(mid(h)*(1:Na)))) ])
h = a0[0]

u = X + (Cd*a0).(X) 
v = a0.(X)

 plot(u,v, xlabel="ξ(x)", ylabel="η(x)", lw=3, guidefontsize=16, tickfontsize=12, titlefontsize=18, legend=false)
 savefig("wave_profile.pdf")