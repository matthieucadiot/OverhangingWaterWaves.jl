using RadiiPolynomial, IntervalArithmetic, LinearAlgebra, Plots, JLD2
using Logging ; disable_logging(Logging.Info)
include("list_fct_branch.jl")

#### Combined existence + geometry analysis of the branch, segment by segment.
#### For each segment we
####   (1) prove the existence of the branch (proof_branch) ;
####   (2) classify the segment :  graph (ξx>0) / overhang (∀m ∃x ξx<0)  and  injective (ξ>-π) / self-intersecting (∀m ∃x ξ<-π) ;
####   (3) if the segment is not cleanly one of these — i.e. a transition occurs inside it — compute a rigorous
####       enclosure Vm of ∂mŨ(m) on the grid together with its error bound ε (‖Vm - ∂mŨ‖_X ≤ ε), save Vm, and
####       display ε so it can be recorded.
####   ξx = D_y a  (eq 10.3) ,  ξ = x + Hh a , (Hh a)_n = coth(nh) aₙ  (eq 10.4).
#### The branch data U_fft_grid_*.jld2 is already saved (in the folder `data` below).
#### The branch can be computed using the following code : save each piece of the branch in extended precision. One file  U_fft_grid_m_<mi0>_<mi>.jld2  is written per segment in the current folder#######  
#save_branch(m0, m1, dm, N, U0, K0, γ, g, h, k ; digits = 128, tol = 1e-30)

K0 = 100 ; N = 12 ; γ = -I"5" ; g = I"1" ; h = I"2"
m0 = -0.1 ; m9 = -1.0 ; dm = 0.1
data = "."                        #### folder containing the branch data (same folder as this script)

dx_grid = 0.002 ; dτ_grid = 0.002           #### x-cell width on (-π,0) ; finest m-cell (the m-direction is adaptive)
nx = ceil(Int, π/dx_grid) ; xnodes = range(-π, 0, length = nx+1)

#### operator giving ξ = x + Hh a ,  (Hh a)_n = coth(nh) aₙ   (eq 10.4)
Hh = LinearOperator(CosFourier(K0,interval(1)), SinFourier(K0,interval(1)),
        [interval.(zeros(K0)) Matrix(Diagonal(coth.(h*interval.(1:K0))))])

#### classify one τ-box (a range of m) over the fixed x-grid : ξx>0 for all x (graph) ? ∃x ξx<0 (overhang) ?
#### ξ>-π for all x (injective, via the trough slope) ? ∃ interior x with ξ<-π (self-intersecting) ?
function classify_box(ξx_cheb, ξmx_cheb, τI, βx, nx, xnodes)
    box_graph = true ; box_overh = false ; box_inj = true ; box_nonphys = false ; incr = true
    for i = 1:nx
        xlo = (i==1) ? inf(-interval(π)) : xnodes[i] ; xI = interval(xlo, xnodes[i+1])
        ξx_eval = ξx_cheb(τI,xI)       + interval(-1,1)*βx
        ξ_eval  = xI + ξmx_cheb(τI,xI) + interval(-1,1)*βx
        ξx_pos = inf(ξx_eval) > 0 ; incr = incr && ξx_pos
        ξx_pos                                     || (box_graph   = false)
        (sup(ξx_eval) < 0)                        && (box_overh   = true)
        (incr || inf(ξ_eval) > sup(-interval(π))) || (box_inj     = false)
        (i>1 && sup(ξ_eval) < inf(-interval(π)))  && (box_nonphys = true)
    end
    return box_graph, box_overh, box_inj, box_nonphys
end

#### NON-RIGOROUS Float64 pre-check (NO interval arithmetic) : does a transition occur in this segment ?
#### We evaluate min_x ξx and the interior min_x ξ at every m-node of the candidate and look for a sign change
#### across the segment. This only DECIDES whether to ask proof_branch for the ∂m (second-derivative) enclosure ;
#### the geometry is proved rigorously a posteriori below, which is what actually certifies the transition.
function segment_has_transition(U_fft, K0, hf, npts_)
    xs  = range(-π, 0, length = 400)
    cth = [coth(k*hf) for k = 1:K0]
    seen_graph = false ; seen_overh = false ; seen_inj = false ; seen_selfint = false
    for n = 1:2*(npts_-1)
        a = [mid(component(U_fft[n],2)[k]) for k = 0:K0]                 #### Float64 cosine coefficients a_0..a_K0
        mξx = Inf ; mξ_int = Inf
        for x in xs
            ξx = 1.0 ; ξ = x
            for k = 1:K0
                ck = cth[k]*a[k+1] ; ξx += 2*k*ck*cos(k*x) ; ξ += 2*ck*sin(k*x)
            end
            mξx = min(mξx, ξx)
            (-π+0.02 < x < -0.02) && (mξ_int = min(mξ_int, ξ))          #### interior min (ξ=-π at the trough is not a self-intersection)
        end
        mξx    > 0  ? (seen_graph = true) : (seen_overh   = true)
        mξ_int > -π ? (seen_inj   = true) : (seen_selfint = true)
    end
    return (seen_graph && seen_overh) || (seen_inj && seen_selfint)     #### both signs present ⇒ a transition is inside the segment
end




for mi = m0-dm : -dm : m9
    mi0 = mi + dm
    fname = joinpath(data, "U_fft_grid_m_$(round(mi0,digits=6))_$(round(mi,digits=6)).jld2")
    U_fft_big = load(fname, "U_fft_grid")
    npts_ = length(U_fft_big) ÷ 2 + 1
    U_fft = [interval.(Float64.(inf.(U_fft_big[n]),RoundDown), Float64.(sup.(U_fft_big[n]),RoundUp)) for n = 1:2*(npts_-1)]

    #### non-rigorous pre-check : does a transition occur in this segment ? (this decides ∂m_control below)
    do_dm = segment_has_transition(U_fft, K0, mid(h), npts_)

    #### (1) existence proof on the segment — run ONCE — asking for the ∂m (second-derivative) enclosure
    #### ONLY on the segments the pre-check flagged as containing a transition.
    println("\n" * "="^74)
    println("EXISTENCE PROOF for the branch on  m ∈ [$mi , $mi0]" *
            (do_dm ? "   (transition flagged — also computing the ∂mŨ enclosure)" : ""))
    println("  Computing the Newton–Kantorovich bounds (Y, Z₁, Z₂) ...")
    println("-"^74)
    Vm_grid = nothing ; ϵ = nothing
    if do_dm
        bounds, Vm_grid, ϵ = proof_branch(K0, N, mi0, mi, h, γ, g, U_fft, U_fft_big ; ∂m_control=1)
    else
        bounds = proof_branch(K0, N, mi0, mi, h, γ, g, U_fft, U_fft_big ; ∂m_control=0)
    end
    r⁻ = bounds.rm
    println("-"^74)
    println("  ✓ EXISTENCE PROVEN on  m ∈ [$mi , $mi0].")
    println("    Newton–Kantorovich bounds :  Y  = $(sup(bounds.Y))")
    println("                                 Z₁ = $(sup(bounds.Z1))")
    println("                                 Z₂ = $(sup(bounds.Z2))")
    println("    Validated radii           :  r⁻ = $(sup(bounds.rm))   (contraction radius / proven accuracy)")
    println("                                 r⁺ = $(inf(bounds.rp))   (largest radius of local uniqueness)")
    println("    ⇒ a unique steady wave U(m) exists in Bᵣ(Ū(m)) for every r ∈ [r⁻, r⁺] and every m ∈ [$mi , $mi0].")

    #### (2) geometry. Chebyshev (in m) coefficients of ξx = D_y a and of (ξ-x) = Hh a.
    println("  Classifying the geometry (graph / overhang / injectivity) on  m ∈ [$mi , $mi0] ...")
    ξx_grid  = [Dy_fct_1D(component(U_fft[n],2),h) for n = 1:2*(npts_-1)]
    ξmx_grid = [Hh*component(U_fft[n],2)           for n = 1:2*(npts_-1)]
    ξx_cheb  = Vector{Sequence{Chebyshev,Vector{Interval{Float64}}}}(undef, K0+1)
    ξmx_cheb = Vector{Sequence{Chebyshev,Vector{Interval{Float64}}}}(undef, K0)
    for n = 1:K0+1
        ξx_cheb[n] = real.(rifft!(complex.([coefficients(ξx_grid[j])[n] for j = 1:2*(npts_-1)]), Chebyshev(N)))
    end
    for n = 1:K0
        ξmx_cheb[n] = real.(rifft!(complex.([coefficients(ξmx_grid[j])[n] for j = 1:2*(npts_-1)]), Chebyshev(N)))
    end
    ξx_cheb  = Sequence(Chebyshev(N)⊗CosFourier(K0,interval(1)), reduce(vcat, [coefficients(ξx_cheb[n])  for n = 1:K0+1]))
    ξmx_cheb = Sequence(Chebyshev(N)⊗SinFourier(K0,interval(1)), reduce(vcat, [coefficients(ξmx_cheb[n]) for n = 1:K0]))

    #### classification over the segment (τ ∈ [-1,1] ↔ m ∈ [mi,mi0]). Defect on ξx and ξ : coth(h) r⁻.
    βx = coth(h)*r⁻
    #### a segment is named for a pure state only if that state holds at EVERY m ; if it holds only at some m,
    #### the segment carries a transition and is named for it.
    #### ADAPTIVE in m : start with one box for the whole segment and bisect a box only when the rigorous
    #### test is inconclusive, down to a finest width dτ_grid. Coarse where the wave is clearly a graph
    #### (resp. overhanging / non-physical), refined only near a transition.
    graph = true ; any_graph = false ; injective = true ; selfint = true ; any_selfint = false ; n_boxes = 0
    stack = [(-1.0, 1.0)]
    while !isempty(stack)
        (τa, τb) = pop!(stack) ; τI = interval(τa, τb)
        bg, bo, bi, bn = classify_box(ξx_cheb, ξmx_cheb, τI, βx, nx, xnodes)
        if ((bg || bo) && (bi || bn)) || (τb - τa ≤ dτ_grid)          #### box resolved (or at finest size) : record it
            bg ? (any_graph = true)   : (graph = false)               #### graph at some m ? / at every m ?
            bn ? (any_selfint = true) : (selfint = false)             #### self-intersecting at some m ? / at every m ?
            bi || (injective = false)                                 #### injective at every m ?
            n_boxes += 1
        else
            τm = 0.5*(τa+τb) ; push!(stack, (τa,τm)) ; push!(stack, (τm,τb))   #### refine near a transition
        end
    end

    label = graph       ? "GRAPH" :
            any_graph    ? "TRANSITION TO OVERHANG" :        #### a graph part at the top, overhanging below
            selfint      ? "NON-PHYSICAL" :                  #### self-intersecting at every m
            any_selfint  ? "TRANSITION TO NON-PHYSICAL" :    #### injective at the top, self-intersecting below
            injective    ? "OVERHANG (injective)" :          #### overhanging, injective everywhere
                           "OVERHANG"
    println("  → Geometry on m ∈ [$mi , $mi0] :  $label    ($(n_boxes) adaptive m-boxes)")

    #### (3) the transition is now established RIGOROUSLY (a posteriori) in this segment. Save the second-derivative
    #### enclosure Vm computed in step (1) : ∂mξx = ∂m∂x ξ is a mixed SECOND DERIVATIVE of the free surface and ∂mξ
    #### a first derivative ; both are recovered from Vm by applying D_y (resp. Hop).
    if label == "TRANSITION TO OVERHANG" || label == "TRANSITION TO NON-PHYSICAL"
        if do_dm
            gname = joinpath(data, "Vm_grid_m_$(round(mi0,digits=6))_$(round(mi,digits=6)).jld2")
            save(gname, "Vm_grid", Vm_grid, "epsilon", ϵ)
            println("   ∂mŨ(m) saved → $gname")
            println("   Rigorous error bound :  ‖Vm - ∂mŨ‖_X ≤ ε = $(sup(ϵ))    (record this value).")
        else
            println("   ⚠ transition rigorously found but the float pre-check missed it — re-run this segment with ∂m_control=1.")
        end
    end
end
