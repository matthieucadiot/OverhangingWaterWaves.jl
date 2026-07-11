# Constructive proofs of a single overhanging wave, a local branch and a global branch of overhanging periodic steady water waves



Table of contents:


* [Introduction](#introduction)
* [Steady water waves with constant vorticity](#steady-water-waves-with-constant-vorticity)
* [Constructive proof of a single overhanging wave](#constructive-proof-of-a-single-overhanging-wave-folder-single-wave-proof)
* [Proof of the local bifurcating branch](#proof-of-the-local-bifurcating-branch-folder-local_branch)
* [Proof of the global branch and its geometric transitions](#proof-of-the-global-branch-and-its-geometric-transitions-folder-main_branch)
* [Installation and reproducibility](#installation-and-reproducibility)
* [Utilisation and References](#utilisation-and-references)
* [License and Citation](#license-and-citation)
* [Contact](#contact)



# Introduction

This Julia code is a complement to the article

#### [[1]](https://arxiv.org/) : "Global bifurcation and the constructive existence of overhanging periodic steady water waves", Matthieu Cadiot and Susanna V. Haziot [ArXiv Link](https://arxiv.org/)

as it provides the rigorous computations supporting the computer-assisted proofs of the paper. Each proof relies on the Newton–Kantorovich theorem: around a numerical approximation $\overline{U}$, one computes a set of rigorous bounds $Y, Z_1, Z_2$ such that, whenever the associated *radii polynomial* admits a root $r$, a unique exact solution is guaranteed within distance $r$ of $\overline{U}$. All floating-point operations are carried out with rigorous outward rounding using the package [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl). The mathematical objects (Fourier and Chebyshev sequence spaces, operators, ...) are built using the package [RadiiPolynomial](https://github.com/OlivierHnt/RadiiPolynomial.jl).

The repository is organised in three folders, mirroring the three constructive results of [[1]](https://arxiv.org/):

- [`single wave proof`](#constructive-proof-of-a-single-overhanging-wave-folder-single-wave-proof) : existence of a single, strongly overhanging wave at a fixed value of the mass flux, together with the verification of its geometric properties;
- [`local_branch`](#proof-of-the-local-bifurcating-branch-folder-local_branch) : existence of the local branch bifurcating from the flat state, which provides the starting point of the continuation;
- [`main_branch`](#proof-of-the-global-branch-and-its-geometric-transitions-folder-main_branch) : existence of the global branch, continued until it self-intersects, together with the rigorous certification of its geometric transitions.


# Steady water waves with constant vorticity

We consider steady, two-dimensional, periodic gravity water waves with constant vorticity $\gamma$, propagating at wave speed $c$ over a flat bed. Using a conformal map from a fixed strip to the (a priori unknown, possibly overhanging) fluid domain, the free-boundary Euler equations are reduced to a single scalar equation posed on the free surface. Writing the surface elevation as a cosine Fourier series $a$, and denoting by $Q$ the Bernoulli constant and by $m$ the relative mass flux, the profiles are the zeros of the nonlinear map $F(U)=0$ with augmented unknown $U=(Q,a)$:

$$
F(U) =
\begin{pmatrix}
a_0 - h \\
\left(m\,\mathbb{D}_y e_0 + \tfrac{1}{2}\gamma\,\mathbb{D}_y(a^2) - \gamma\,a\,\mathbb{D}_y a\right)^2 - (Q - 2g\,a)\left((\mathbb{D}_x a)^2 + (\mathbb{D}_y a)^2\right)
\end{pmatrix}
= 0.
$$

Here $g$ is the gravitational constant, $h$ the (conformal) mean depth, $e_0$ the constant identity sequence, and $\mathbb{D}_x, \mathbb{D}_y$ are the Fourier multiplier operators induced by the conformal map: $\mathbb{D}_x$ is spatial differentiation (symbol $in$ on mode $n$) and $\mathbb{D}_y$ has symbol $n\coth(nh)$ (and $1/h$ at $n=0$). The first component fixes the conformal depth, while the second is the Bernoulli dynamic boundary condition.

The physical profile is reconstructed through $\xi_x = \mathbb{D}_y a$ and $\xi = x + \mathbb{H} a$, where $(\mathbb{H} a)_n = \coth(nh)\,a_n$. A profile is a **graph** when $\xi_x>0$ everywhere, is **overhanging** when $\xi_x<0$ somewhere, and is **non-physical** (self-intersecting on the trough line) when $\xi<-\pi$ at some interior point. These are precisely the quantities certified in the geometric part of the proofs.

Throughout the paper and the code, the physical parameters are fixed to $\gamma=-5$ (vorticity), $g=1$ (gravity) and $h=2$ (conformal depth); the relative mass flux $m$ is used as the continuation parameter.


# Constructive proof of a single overhanging wave (folder `single wave proof`)

The code `main_proof.jl` provides the computer-assisted proof of existence of a single, strongly overhanging wave at the fixed mass flux $m=-0.85$ (Section 8 of [[1]](https://arxiv.org/)). Starting from a numerical approximation (stored in `U.jld2`, the corresponding mass flux in `m.jld2`), the script proceeds in three steps.

First, it assembles the Jacobian $DF(\overline{U})$ together with an approximate inverse $A$. The construction of $A$ combines a finite-dimensional numerical inverse with an *analytic tail correction*, so that $A$ acts as an accurate inverse on the whole (infinite-dimensional) sequence space and not only on the truncated modes.

Second, it computes rigorously the Newton–Kantorovich bounds $Y$, $Z_1$ and $Z_2$ and evaluates the radii polynomial. When the proof succeeds, the script prints the three bounds and the validated contraction radius $r^-$: a unique exact wave lies within distance $r^-$ of the approximation.

Third, it verifies the geometric and non-degeneracy properties of the proven wave: the surface is positive ($\eta>0$), monotone on a half-period, strictly overhanging (there is a point where $\xi_x<0$), and injective (the profile does not self-intersect). The certified profile is finally plotted over two periods and exported to `wave_profile.pdf`.

All the specialised routines — the operators $\mathbb{D}_x,\mathbb{D}_y$, the weighted norms, and the explicit second-derivative bounds used for the monotonicity argument — are contained in `list_fct_proof.jl`.


# Proof of the local bifurcating branch (folder `local_branch`)

Nontrivial waves bifurcate from the trivial branch of flat states through a pitchfork bifurcation. The purpose of this folder is to prove, constructively and quantitatively, the existence of this **local branch** and to advance far enough along it to provide a rigorous starting point for the global continuation (Section 5 of [[1]](https://arxiv.org/)).

The branch is parameterised by $\epsilon$, with $m = \overline{m} - \epsilon^2$ ($\overline{m}$ being the bifurcation value of the mass flux), and the approximate branch is represented by a Chebyshev polynomial in $\epsilon$ whose coefficients are themselves cosine Fourier series. The script `construction_branch.jl` builds this approximation (Newton iterations on a Chebyshev grid) and saves it to `U_fft.jld2`, together with the end point of the branch (`U_end.jld2`, `m_end.jld2`) that is later used to launch the global continuation.

The main script `main_local_bif.jl` then:

1. loads the approximate branch and runs `proof_branch`, which computes the Newton–Kantorovich bounds $(Y, Z_1, Z_2)$ uniformly along the branch and returns the smallest and largest validated radii $r^-$ and $r^+$. A unique branch of exact solutions is enclosed in the tube of radius $r$ around the approximation for every $r\in[r^-,r^+]$;
2. proves, via `proof_geometry`, that every wave of the local branch is a graph ($\xi_x>0$ on the half-period for all $\epsilon\in[0,\epsilon_0]$) and that the profile is monotone for $\epsilon$ small enough.

The functions used here — and reused by the global-branch folder — are collected in `list_fct_branch.jl`.


# Proof of the global branch and its geometric transitions (folder `main_branch`)

This folder contains the central result of [[1]](https://arxiv.org/): the rigorous continuation of the branch from the local regime ($m=-0.1$) down to $m=-1.0$, and the proof that, along the way, the branch transitions from graphs to overhanging profiles and finally terminates at a self-intersecting (non-physical) wave. This resolves, in the parameter regime under consideration, the topological termination conjecture of Constantin, Strauss and Varvaruca (Sections 6 and 7 of [[1]](https://arxiv.org/)).

The interval $m\in[-1,-0.1]$ is split into nine segments of width $0.1$. For each segment, an approximate solution — again a Chebyshev polynomial in $m$ with cosine-Fourier coefficients — is stored in a file `U_fft_grid_m_<mi0>_<mi>.jld2` (for instance `U_fft_grid_m_-0.1_-0.2.jld2`). Because these nine files together exceed GitHub's file-size limit, they are shipped compressed as `data_branch.zip`, which must be unzipped inside `main_branch/` before running the proof, so that the nine `.jld2` files sit directly next to `main_existence_geometry.jl`. The end-point data `U_end.jld2` and `m_end.jld2` are inherited from the local branch.

The single driver `main_existence_geometry.jl` performs, segment by segment:

**1. Existence.** It calls `proof_branch` to compute the Newton–Kantorovich bounds $(Y, Z_1, Z_2)$ and the validated radii $r^-, r^+$, establishing a unique smooth branch of exact solutions on the segment. Adjacent segments are then *glued* together through the overlap of their uniqueness balls, producing a single continuous branch over the whole interval. For every segment the script displays a header, the three bounds, and the two radii.

**2. Geometric classification.** From the same candidate it forms the Chebyshev–Fourier expansions of $\xi_x=\mathbb{D}_y a$ and $\xi = x + \mathbb{H} a$, and classifies the segment — adaptively refining in $m$ only where needed — as `GRAPH`, `OVERHANG`, injective or `NON-PHYSICAL`, always adding the rigorous a-posteriori defect $\coth(h)\,r^-$ to the candidate before taking the sign.

**3. Certifying the transitions.** On the (single) segment where the geometry changes, the script additionally asks `proof_branch` for a rigorous enclosure of the parameter derivative $\partial_m \tilde U(m)$ together with its error bound $\epsilon$. From this enclosure one controls the mixed derivative $\partial_m \xi_x = \partial_m\partial_x \xi$, which certifies that $\min_x \xi_x$ is *strictly monotone* in $m$ across the gap. Consequently the graph $\to$ overhang transition — and, by the same argument applied to $\xi$ versus $-\pi$, the injective $\to$ self-intersecting transition — occurs **exactly once**. The corresponding enclosure is saved as `Vm_grid_m_<...>.jld2`.

All the operators, the computation of the bounds, and the $\partial_m$-control machinery are contained in `list_fct_branch.jl`.


# Installation and reproducibility

The exact versions of every package used to run the proofs are pinned in the `Project.toml` and `Manifest.toml` files at the root of the repository. The proofs were run with **Julia 1.12** and, in particular, `RadiiPolynomial` v0.9.11, `IntervalArithmetic` v1.0.10, `JLD2` v0.6.3 and `Plots` v1.41.3. To recreate the environment, start Julia from the repository root and instantiate it:

```julia
using Pkg
Pkg.activate(".")     # use the environment shipped with this repository
Pkg.instantiate()     # install the exact versions recorded in Manifest.toml
```

or, equivalently, from a terminal:

```
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

The global-branch data is shipped compressed, to stay within GitHub's file-size limit; unzip it inside `main_branch/` before running that proof:

```
cd main_branch && unzip data_branch.zip     # restores the nine U_fft_grid_m_*.jld2 files
```

Because the rigorous bounds depend on the exact behaviour of `IntervalArithmetic` (and of `RadiiPolynomial`, which pins it), instantiating the shipped `Manifest.toml` is what guarantees an identical reproduction of the computer-assisted proofs.


# Utilisation and References

To reproduce a given result, download the content of the corresponding folder into a single directory and run its main script from that directory:

- `single wave proof/main_proof.jl` — single overhanging wave and its geometry;
- `local_branch/main_local_bif.jl` — local bifurcating branch (`construction_branch.jl` regenerates the approximation, if needed);
- `main_branch/main_existence_geometry.jl` — global branch and its geometric transitions (first run `unzip data_branch.zip` in that folder to restore the nine `U_fft_grid_m_*.jld2` data files).

The files `list_fct_proof.jl` and `list_fct_branch.jl` contain all the functions needed to build the operators and compute the rigorous bounds, while the `.jld2` files hold the pre-computed numerical approximations. The codes are written for the parameter values of [[1]](https://arxiv.org/) ($\gamma=-5$, $g=1$, $h=2$), but the routines are general and compute the required bounds for other admissible values as well.

The codes are built using the following packages:
- [RadiiPolynomial](https://github.com/OlivierHnt/RadiiPolynomial.jl)
- [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl)
- [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/)
- [JLD2](https://github.com/JuliaIO/JLD2.jl)
- [Plots](https://github.com/JuliaPlots/Plots.jl) (only used by `main_proof.jl` to export the wave profile)


# License and Citation

This code is available as open source under the terms of the [MIT License](http://opensource.org/licenses/MIT).

If you wish to use this code in your publication, research, teaching, or other activities, please cite it using the following BibTeX template:

```
@software{OverhangingWaterWaves.jl,
  author = {Matthieu Cadiot},
  title  = {OverhangingWaterWaves.jl},
  url    = {https://github.com/matthieucadiot/OverhangingWaterWaves.jl},
  note = {\url{https://github.com/matthieucadiot/OverhangingWaterWaves.jl}},
  year   = {2026}
}
```

# Contact

You can contact me at :

matthieu.cadiot@polytechnique.edu
