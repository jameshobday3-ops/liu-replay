# Window 1.5 — monotone block-Krylov inverse hierarchy

**Status.** Exact operator lemma for the auxiliary inverse-majorant construction. It does not prove a spectral gap at (L=1.5) and does not prove RH. It explains why enlarging the inverse trial space is a principled continuation of the certified comparison method.

## 1. Setup

Let (M) be bounded self-adjoint on a Hilbert space (H), with
[
0<mIpreceq Mpreceq bI,
]
and let (E:\mathbb C^N\to H) be an isometry onto the finite polynomial space used in the comparison. Put
[
G=E^*M^{-1}E.
]

For any trial map (Y:\mathbb C^N\to H), define its residual
[
R_Y=E-MY
]
and the standard residual majorant
[
\mathcal W(Y)
=
E^*Y+Y^*E-Y^*MY+m^{-1}R_Y^*R_Y.
]

A direct completion of squares gives
[
\boxed{
\mathcal W(Y)-G
=
R_Y^*(m^{-1}I-M^{-1})R_Y
\succeq0.
}
\tag{1}
]

Thus every such trial map gives an upper bound on the compressed inverse.

## 2. Monotonicity under nested trial spaces

Let
[
\mathcal V_0\subseteq\mathcal V_1\subseteq\cdots\subset H
]
be nested finite-dimensional trial spaces, and let (Y_r) minimize the quadratic majorant over maps with
[
\operatorname{ran}(Y)\subseteq\mathcal V_r.
]

Then
[
\boxed{
G\preceq \mathcal W(Y_{r+1})\preceq \mathcal W(Y_r).
}
\tag{2}
]

The second inequality is simply the variational inclusion: every map admissible at rank (r) remains admissible at rank (r+1). Equivalently, the block Schur update obtained when a new residual direction is appended is positive semidefinite in the inverse target and therefore improves the lower comparison in the useful direction.

This is the abstract form of the rank-one and rank-two residual enrichments already used at (L=1.5).

## 3. First block-Krylov space

The natural first enlargement is
[
\mathcal V_1
=
\operatorname{ran}(E)+\operatorname{ran}(ME).
]

For a seed (a\in\mathbb C^N), decompose
[
MEa=Eb+v,qquad b=\Gamma^{-1}E^*MEa,qquad v=(I-\Pi)MEa.
]

Appending (v) therefore adds a direction in the first block-Krylov complement
[
(I-\Pi)ME.
]

The adaptive residual construction used in the Window-15 computation selects seed vectors (a_j) from surviving negative comparison directions and adds
[
v_j=(I-\Pi)MEa_j.
]
After (Gamma)-orthogonalizing the (a_j), this produces a nested sequence inside the same first block-Krylov space. Its acceptance still uses fresh directed (M^3) columns, (M^4) Gram entries and a full Arb inertia calculation; the floating/high-precision eigendirection used to propose (a_j) has no evidential role.

## 4. Complete Krylov hierarchy

More generally let
[
\mathcal K_k(M,E)
=
\operatorname{span}\{E,ME,M^2E,\ldots,M^kE\}.
]

For a polynomial (p_k) of degree at most (k), the trial
[
Y=p_k(M)E
]
has residual
[
R_Y=q_{k+1}(M)E,qquad q_{k+1}(x)=1-xp_k(x),qquad q_{k+1}(0)=1.
]

Consequently, if
[
\epsilon_k=\max_{x\in[m,b]}|q_{k+1}(x)|,
]
then from (1)
[
0\preceq \mathcal W(Y)-G
\preceq
(m^{-1}-b^{-1})\epsilon_k^2 I.
\tag{3}
]

(The coarser constant (m^{-1}epsilon_k^2) is sufficient if the upper spectral bound is not used.)

The scaled Chebyshev residual polynomial gives
[
\epsilon_k
\le
\frac{2\rho^{k+1}}{1+\rho^{2k+2}}
<2\rho^{k+1},
qquad
\rho=
\frac{\sqrt\kappa-1}{\sqrt\kappa+1},
qquad
\kappa=\frac bm.
\tag{4}
]

Hence
[
\boxed{
\|\mathcal W_k-G\|
\le
4(m^{-1}-b^{-1})\rho^{2k+2},
}
\tag{5}
]
for the optimal majorant over the full block-Krylov space, since it can only improve on the explicit Chebyshev trial.

## 5. Window-15 constants

For the certified adaptive source at (L=1.5),
[
\frac12 I\preceq M\preceq9.1I.
]
Thus
[
\kappa=18.2,
qquad
\rho=
\frac{\sqrt{18.2}-1}{\sqrt{18.2}+1}
\approx0.6202156044720571.
]

Therefore the inverse-majorant hierarchy converges geometrically to the exact compressed inverse. This does **not** by itself say how many Krylov layers are required to close the physical parity gap: the final comparison may be extremely close to singular, and the downstream (K)-congruence can amplify inverse error. It does prove that failure of a low-rank residual enrichment is a limitation of that finite trial space, not an intrinsic obstruction of the inverse-majorant principle.

## 6. Current computational use

The present adaptive computation deliberately exhausts a small number of directions inside the first block
[
(I-\Pi)ME
]
before paying for the next Krylov power. Each new seed requires only targeted (M^3) and (M^4) data.

If this first block remains insufficient, the next principled step is to add directions from
[
(I-\Pi)M(I-\Pi)ME,
]
which enters the (k=2) block-Krylov layer and requires the corresponding higher moments. The hierarchy (2)–(5) then continues unchanged.

## 7. Logical role in the RH programme

The fixed-window hierarchy is a calibration device for the more important Feshbach target. Once a simple isolated even ground state is obtained, the same inverse-weighted philosophy is applied to the prolate candidate (v_L):
[
\eta_L=
\|(C_L-U_L)^{-1}Q_LA_Lv_L\|.
]

An all-window proof still requires a scale-uniform coercivity/shorted-resolvent theorem and quantitative decay of (eta_L). Nothing in this note extrapolates a finite (L=1.5) calculation to unbounded windows.
