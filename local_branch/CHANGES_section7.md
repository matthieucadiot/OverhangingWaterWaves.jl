# Update of `list_fct_branch.jl` to the new Section 7 bounds

All changes are confined to `list_fct_branch.jl` (functions `G_grid` and `proof_branch`).
`main_local_bif.jl` is **unchanged**: `proof_branch`'s signature and call are identical.

Each edit is marked in the code with a `#### NEW (...)` comment. Mapping to the PDF:

## 1. Approximate inverse `A_K` — Sec. 7.2, eq. (7.14) / p.29
- The tail operator `A∞ = π^{>K1} D_y^{-1} M_w` (`Dy_Mw_grid`) is now built **before** `A_K`,
  and for **all** grid points (including the bifurcation index `n0`).
- `A_K` is now the numerical inverse of
  `M = π^{≤K1} DG (I − A∞ DG) π^{≤K1}`  (variable `M_grid0`),
  instead of the plain Galerkin inverse `inv(π DG π)`.
  Only the a-row of `M` gets the `A∞`-correction (the Q-row reads mode 0, which `A∞ DG` does not feed).
- The Chebyshev/DFT re-expansion of `A_K` is unchanged.

## 2. `Z1` block — Lemma 7.3, eq. (7.16)
- **Z10** (new combined form): `‖I − A_K M‖` with the `M` above
  ( = `‖π^{≤K1} − A_K P DF (π^{≤K1} − A∞ DF) P π^{≤K1}‖` ).
- **Z11**: unchanged.
- **Z12**: this is the **former `Z14`** = `‖A_K DF π^{>K1}‖`. The old `Z12`, `Z13` are removed.
- **Z∞**: first-term coefficient changed `(K1+2)/(K1+1)^2  →  1/(K1+1)`. Second term unchanged.
- **δ**: decay rate changed `e^{−2 d K0}  →  e^{−2 d K1}` (both occurrences). `d = h`.
- **Z1** = `max{ Z10 + Z11 , (1 + Z12) Z∞ }`   (was `max{Z10+Z11+Z12, (1+Z14)Z∞+Z13}`).

## 3. `Y` block — Lemma 7.4 (updated to Water_Waves_2D-25.pdf)
- `G_grid` now also returns the raw residual `f = F(V, m̄−ε²)` (`F_raw`; `f = 0` at `n0`).
- **Term 1** (`AGG`): `(1/ε²)‖A_K (I − DG·A∞) P F‖_H1 = ‖A_K (I − DG·A∞) G‖_H1`.
  Built exactly like the `M_grid`/`A_K` construction: only the a-row gets the `DG·A∞` correction
  (`A∞ G` is a tail, so the Q-row reads mode 0 → 0).
- **Term 2**: `(K1+2)/(K1+1) ( ‖π^{k>K1}(w*f)‖_ℓ1 + ‖π^{k<−K1}(w̄*f)‖_ℓ1 )`,
  the one-sided tails of the convolution `w*f` (`w̄ = conj(w) = reverse(w)` since `w` has real
  coefficients). No `Z12` appears in the new formula.
- Note: the previous draft (`-26.pdf`) had the older Y (`‖A_K G‖ + (K1+2)/(K1+1)(‖w‖·Z12·‖π^{>K1}f‖ + ‖w*f‖)`);
  this is now replaced by the `-25.pdf` version above.

## 4. `Z2` block — Lemma 7.5
- New term **Z_{D²,1}** (`ZD21`) = `sup ‖A P D²F(V,m̄−ε²)(Z̄,Z̄)‖_H1`,
  computed with the un-rescaled `A_K` (the tail `A∞` vanishes on these low modes).
- **Z_A** (`norm_A`): unchanged formula, with `Z14 → Z12`.
- **Z_{D²,2}** (`ZD2`): factor **2** added in front of `‖Q̄ e0 − 2g ā‖` (third term).
- **Z_{D³}** (`ZD3`): now the bracket `Z_{D³}(r0)` **without** the leading `r0` (applied below).
- **Z2** = `max{ Z_{D²,1}, Z_A·Z_{D²,2}, ε0·Z_A·Z_{D²,2} } + Z_A·Z_{D³}(r0)·r0`
  (was `Z2 = norm_A·(ZD2 + ZD3)`).

## 5. Final radii-polynomial check
- Unchanged (uses `Z1, Z2, Y, r0`).

## Notes for verification
- `Z̄` is encoded as `Zb4` with `component(Zb4,2)[1]=1` (fundamental cos mode), matching the
  existing `Z̄` in `Ḡ`/`DḠ`.
- The new code reuses the existing helpers (`D2F`, `Dy_op_inv`, `D0_op`, `norm_op_grid`,
  `norm_seq_grid`, `cos2exp`, …) and the same interval/Chebyshev machinery as before.
- Not run here (no Julia in this environment): bracket and block-structure balance were checked
  programmatically, but please run the proof once to confirm numerical execution.
