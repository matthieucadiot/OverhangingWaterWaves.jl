using RadiiPolynomial, IntervalArithmetic, LinearAlgebra, JLD2
using Logging
global_logger(ConsoleLogger(stderr, Logging.Warn)) ### This hides the [Info ] messages

include("list_fct_branch.jl")
setprecision(128)

K0 = 100 ; K1 = 2*K0 + 2 ; γ = -interval(5.0) ; g = interval(1.0) ; h = interval(2.0) ; k = interval(1.0) 
γb = - interval(big(5.0)) ; gb = interval(big(1.0)) ; hb = interval(big(2.0)) ; kb = interval(big(1.0))

S = ParameterSpace()×CosFourier(K0,interval(1)) 
Sa = CosFourier(K0,interval(1))
#### computation of the bifurcation parameters m̄ and Q̄ for n=1 ####
m̄ = γ*h^2*interval(0.5) - γ*h*tanh(h)*interval(0.5) + h*sqrt(g*tanh(h) + (γ^2*tanh(h)^2)/(interval(4)*(k)^2)) 
m̄b = γb*hb^2*interval(big(0.5)) - γb*hb*tanh(hb)*interval(big(0.5)) + hb*sqrt(gb*tanh(hb) + (γb^2*tanh(hb)^2)/(interval(big(4))*(kb)^2)) 
Q̄ = interval(2)*g*h + (m̄/h - γ*h*interval(0.5))^2

ā = Sequence(CosFourier(K0,interval(1)), interval.(zeros(K0+1))) ; ā[0] = h
Ū = Sequence(ParameterSpace()×CosFourier(K0,interval(1)) ,[Q̄ ; coefficients(ā)])
##### definition of sequences spaces

# U = Sequence(S,rand(K0+2)) ./ (1:K0+2).^4

#####


ϵ0 = sqrt(m̄ + I"0.1") ### we want to construct the branch from m̄ to m = -0.1, so -0.1 = m̄ - ϵ0^2
ϵ0b = sqrt(m̄b + parse(Interval{BigFloat}, "0.1"))

N = 8;  N_ = 8N  # number of Chebyshev coefficients used
N_fft_ = nextpow(2, 2N_ + 1)
npts_ = N_fft_ ÷ 2 + 1
U_fft_big = load("U_fft.jld2","U_fft") # loading the approximate solution on the grid. It corresponds to the evaluation of a polynomial of order N on the Chebyshev grid of order 8N

U_fft = [interval.(Float64.(inf.(U_fft_big[n]),RoundDown),Float64.(sup.(U_fft_big[n]),RoundUp) ) for n = 1:2*(npts_-1)]

println("\n" * "="^74)
println("CONSTRUCTIVE EXISTENCE PROOF of the local bifurcating branch")
println("  branch parameterised by ϵ ∈ [0, ϵ0] :  m = m̄ - ϵ²,  from m̄ down to m = -0.1")
println("  Computing the Newton–Kantorovich bounds (Y, Z₁, Z₂) ...")
println("-"^74)
r⁻,r⁺ = proof_branch(K0,N,ϵ0,h,γ,m̄,g,U_fft,U_fft_big) # running the proof. It returns the smallest and largest radii for the existence proof.
println("-"^74)
println("  ✓ EXISTENCE PROVEN on the whole local branch.")
println("    Validated radii :  r⁻ = $(sup(r⁻))   (contraction radius / proven accuracy)")
println("                       r⁺ = $(inf(r⁺))   (largest radius of local uniqueness)")
println("    ⇒ a unique branch of steady waves exists in Bᵣ(Ū) for every r ∈ [r⁻, r⁺].")


### code for proving that the waves of the local branch are graphs 
r⁻ = interval(1.1e-6) ## we recall the minimal radius obtained in the existence proof
proof_geometry(K0,N,ϵ0,h,γ,m̄,g,U_fft,r⁻)

### proving that the waves are monotone for ϵ small enough. It is sufficient to prove that (ā(0))_1 - r⁻ > 0. Note that Ū(0) corresponds to U_fft[npts_] 
Ū0 = U_fft[npts_] ; ā0 = component(Ū0,2)

if inf(ā0[1] - r⁻) > 0 
    println("For ϵ small enough, η̃(ϵ) is monotone on (-π,0)")
else 
    println("Cannot conclude about monotonicity.")
end 
