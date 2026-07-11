function cos2exp(a)
    N = size(a)[1]-1
    return Sequence(Fourier(N,interval(1)), vec([reverse(a[1:N]); a[0:N]]))
end

function sin2exp(a)
    N = size(a)[1]
    return Sequence(Fourier(N,interval(1)), vec([im*reverse(a[1:N]); interval(0) ; -im*a[1:N]]))
end


function Dy_fct_1D(U,d::Interval{Float64})
    Sc = space(U)
    K = order(Sc)
    dn = interval.(1:K)
    dn = [interval(1)/d ; dn.*coth.(d*dn)]
    return dn.*U
end 



function Dy_fct_1D(U,d::Float64)
    Sc = space(U)
    K = order(Sc)
    dn = (1:K)
    dn = [(1)/d ; dn.*coth.(d*dn)]
    return dn.*U
end 



# function D0_fct(U::Sequence{CartesianProduct{Tuple{ParameterSpace, CosFourier{Interval{Float64}}}}, Vector{Interval{Float64}}})
#     Sc = space(U) ; K = order(Sc[2])
#     dn = [interval(1) ; interval(1) ; interval(2)*interval.(2:K+1)]
#     return dn.*U
# end 


function D0_fct(U)
    Sc = space(U) ; K = order(Sc)
    dn = [interval(1) ; interval(2)*interval.(2:K+1)]
    return dn.*U
end 



#### Rigorous upper bound for the error on the second derivative of the solution :
####     ‖ D_x² ã - ā_xx ‖_X ≤ bound_error_axx(...)         (Lemma 7.12 + Remark 7.13)
#### a,Q : candidate (CosFourier) and parameter ; w,δ : Lemma 7.6 ; r0 : validation radius.
function bound_error_axx(a, Q, m, h, γ, g, w, δ, r0)
    cth = coth(h)
    e0 = Sequence(space(a), interval.(zeros(size(a)[1]))) ; e0[0] = interval(1)

    #### T₀(ā)  (Lemma 7.12) : even part (cosine) + odd part (sine), then exp representation
    Dxa = Derivative(1)*a ; Dya = Dy_fct_1D(a,h) ; Dya2 = Dy_fct_1D(a^2,h)
    P1 = interval(0.5)*γ*Dya2 - γ*a*Dya                          # cosine
    P2 = m*Dy_fct_1D(e0,h) + interval(0.5)*γ*(Dya2 - a*Dya)      # cosine
    P3 = g*Dxa*(Dxa^2 + Dya^2)                                  # sine
    qa = Q*e0 - interval(2)*g*a
    T0 = cos2exp(qa*P1*P2) - sin2exp(qa*P3)
    w2 = w*reverse(w) ; Tbar = w2*T0 ; nT = norm(Tbar,1)        # |w|², T̄(ā)=|w|²T₀(ā), ‖T̄(ā)‖_{ℓ1}

    aX = norm(D0_fct(a)) ; al1 = norm(a,1)                    # ‖ā‖_X , ‖ā‖_{ℓ1}

    #### ε₁  (Lemma 7.12)
    ε1 = δ + norm(w,1)*( max(interval(1),interval(2)*g)*(interval(1)+cth)*aX
           + (abs(Q)+r0+interval(2)*g*(al1+r0))*(interval(1)+cth) )*r0

    #### ‖T₀(ã)-T₀(ā)‖_{ℓ1} via the telescoping product bound of Lemma 18.3 (Remark 7.13).
    #### For each factor X(a) of T₀ we carry the pair (β,d):
    ####     β  = an upper bound on ‖X(ã)‖_{ℓ1}        (the size of the factor at the true solution)
    ####     d  = an upper bound on ‖X(ã)-X(ā)‖_{ℓ1}   (the "error term" of that factor)
    #### β is always  ‖X(ā)‖_{ℓ1} (computable) + d ; the error d propagates through products by
    #### Lemma 18.3 :  d_{u*v} = d_u β_v + β_u d_v .  We use ‖ã-ā‖_X ≤ r0 and the operator bounds
    #### ‖D_x b‖_{ℓ1} ≤ ‖b‖_X  and  ‖D_y b‖_{ℓ1} ≤ coth(h)‖b‖_X.
    dX2     = interval(2)*(aX+r0)*r0                            # ‖ã²-ā²‖_X : Lemma 18.3 on a*a in X, 2 δ β = 2 r0 (‖ā‖_X+r0)
    βa,da   = al1+r0, r0                                        # factor  a   : d = ‖ã-ā‖_{ℓ1} ≤ r0
    βx,dx   = norm(Dxa,1)+r0,        r0                         # factor  D_x a : d = ‖D_x(ã-ā)‖_{ℓ1} ≤ ‖ã-ā‖_X ≤ r0
    βy,dy   = norm(Dya,1)+cth*r0,    cth*r0                     # factor  D_y a : d ≤ coth(h)‖ã-ā‖_X ≤ coth(h) r0
    βy2,dy2 = norm(Dya2,1)+cth*dX2,  cth*dX2                    # factor  D_y(a²) : d ≤ coth(h)‖ã²-ā²‖_X = coth(h) dX2
    βay,day = βa*βy,                 da*βy+βa*dy                # factor  a*(D_y a)        (Lemma 18.3, two factors a and D_y a)
    βP1,dP1 = abs(γ)*(interval(0.5)*βy2+βay), abs(γ)*(interval(0.5)*dy2+day)   # P1 = γ/2 D_y(a²) - γ a*(D_y a)  (sum of the two factors above)
    βP2,dP2 = abs(m)/h+abs(γ)*interval(0.5)*(βy2+βay), abs(γ)*interval(0.5)*(dy2+day)  # P2 = m/h + γ/2 (D_y(a²) - a*(D_y a)) ; the constant m/h has zero error
    β12,d12 = βP1*βP2,               dP1*βP2+βP1*dP2            # P1*P2                    (Lemma 18.3, factors P1 and P2)
    βsq,dsq = βx^2+βy^2,             interval(2)*(βx*dx+βy*dy)  # (D_x a)² + (D_y a)²      (each square: d = 2 β d)
    βP3,dP3 = g*βx*βsq,              g*(dx*βsq+βx*dsq)          # P3 = g (D_x a)*((D_x a)²+(D_y a)²)  (factors D_x a and the sum of squares)
    βbr,dbr = β12+βP3,               d12+dP3                    # bracket = P1*P2 - P3     (triangle inequality on the difference)
    βqa,dqa = abs(Q)+interval(2)*g*al1+(interval(1)+interval(2)*g)*r0, (interval(1)+interval(2)*g)*r0   # factor  Q-2ga : d ≤ |Q̃-Q̄| + 2g‖ã-ā‖ ≤ (1+2g) r0
    dT0     = dqa*βbr + βqa*dbr                                 # T₀ = (Q-2ga)*bracket : only the error term d_{T₀} is needed  (Lemma 18.3)

    #### ε₂  (Lemma 7.12 + Remark 7.13)
    ε2 = (norm(w2,1)*dT0 + (interval(2)*ε1+ε1^2)*nT) / (interval(1)-interval(2)*ε1-ε1^2)

    #### bound (7.47)
    return aX*(interval(1)+cth)*ε2 + nT*(interval(1)+cth)*r0
end




#### T₀(ā)  (Lemma 7.12), returned in the exponential (full Fourier) representation.
function T0(a, Q, m, h, γ, g)
    e0 = Sequence(space(a), interval.(zeros(size(a)[1]))) ; e0[0] = interval(1)
    Dxa = Derivative(1)*a ; Dya = Dy_fct_1D(a,h) ; Dya2 = Dy_fct_1D(a^2,h)
    P1 = interval(0.5)*γ*Dya2 - γ*a*Dya                          # cosine
    P2 = m*Dy_fct_1D(e0,h) + interval(0.5)*γ*(Dya2 - a*Dya)      # cosine
    P3 = g*Dxa*(Dxa^2 + Dya^2)                                  # sine
    qa = Q*e0 - interval(2)*g*a
    return cos2exp(qa*P1*P2) - sin2exp(qa*P3)
end

#### T̄(ā) = |w|² * T₀(ā)   (Lemma 7.12).   |w|² = w * w⋆ = w * reverse(w)  (w has real coefficients).
Tbar(a, Q, m, h, γ, g, w) = (w*reverse(w)) * T0(a, Q, m, h, γ, g)

#### periodic Hilbert transform on the strip  H_h : (H_h v)_n = -i tanh(n h) v_n   (Section 12).
Hh(v, h) = Sequence(space(v), [ -im*tanh(interval(k)*h) for k = -order(v):order(v) ] .* coefficients(v))

#### approximate second derivative  ā_xx = T̄(ā)*(D_x ā) + (D_y ā)*(H_h T̄(ā))   (Lemma 7.12),
#### returned in the exponential representation.
function axx(a, Q, m, h, γ, g, w)
    Tb = Tbar(a, Q, m, h, γ, g, w)
    v  = Tb*sin2exp(Derivative(1)*a) + cos2exp(Dy_fct_1D(a,h))*Hh(Tb, h)                # exponential representation
    K  = order(v)
    return Sequence(CosFourier(K,interval(1)), real.(v[0:K]))  # even real part → cosine series
end
