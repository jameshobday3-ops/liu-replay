# RH shorted-resolvent theorem: one-vector certificate for the two missing CCM steps

**2 November 2026 · James Hobday research programme**

## Purpose

This note isolates a sufficient operator theorem that combines the two missing
steps identified by Connes–Consani–Moscovici for their localized Weil route:

1. the lowest eigenvalue is simple and its eigenfunction is even;
2. the prolate candidate is sufficiently close to that ground eigenfunction.

The theorem does **not** assume that the lowest eigenvalue is positive.

## 1. Setting

Let (A) be a lower-bounded self-adjoint operator with compact resolvent on a
Hilbert space (mathcal H). Assume that (A) commutes with an involutive
parity operator (R).

Let (vinoperatorname{Dom}(A)) be a unit even vector,
(Rv=v), and define
[
P=|vanglelangle v|,qquad Q=I-P.
]
Write
[
mu=langle v,Avangle,qquad b=QAv.
]
Let (C) be the self-adjoint operator on (Qmathcal H) associated with the
restriction of the closed quadratic form of (A) to (Qmathcal H).

Assume that for real numbers (U<c),
[
mule U<c,qquad Cge cI
]
in the quadratic-form sense.

## 2. Theorem

Under the assumptions above:

1. (A) has exactly one eigenvalue below (c), counted with multiplicity.
   Denote it by (lambda_1). It is therefore simple and isolated.

2. A normalized eigenvector (g) for (lambda_1) is even.

3. One has the exact Feshbach identity
[
mu-lambda_1
=
leftlangle b,(C-lambda_1)^{-1}bightangle.
]

4. Choosing the phase of (g) so that (langle v,gangle>0),
[
oxed{
|g-v|
le
|(C-lambda_1)^{-1}b|
le
|(C-U)^{-1}b|.
}
]

Thus the single quantity
[
oxed{
eta:=|(C-U)^{-1}QAv|
}
]
controls the candidate-to-ground-state distance while the complement bound
(Cge cI>UI) simultaneously proves the simple-even property.

## 3. Proof

Since (v) is a unit trial,
[
lambda_1(A)le mule U<c.
]

For every two-dimensional subspace (Esubsetmathcal H), there is a nonzero
(fin Ecap v^perp). Hence
[
sup_{0
e hin E}rac{q_A[h]}{|h|^2}
ge
rac{q_A[f]}{|f|^2}
ge c.
]
The min–max principle therefore gives
[
lambda_2(A)ge c.
]
Since (lambda_1<clelambda_2), the eigenvalue below (c) is unique,
simple and isolated.

Because (A) commutes with parity and (lambda_1) is simple, its
eigenvector has definite parity. Every odd vector lies in (v^perp), since
(v) is even. The complement bound therefore places every odd Rayleigh
quotient at or above (c). The eigenvector at (lambda_1<c) cannot be odd,
so it is even.

Write the normalized eigenvector as
[
g=alpha v+h,qquad hin Qmathcal H.
]
If (alpha=0), then (gin Qmathcal H) and its Rayleigh quotient is at
least (c), contrary to (lambda_1<c). Thus (alpha
e0). The (Q)-part
of the eigenvalue equation gives
[
(C-lambda_1)h=-alpha b,
]
so
[
h=-alpha(C-lambda_1)^{-1}b.
]
The (P)-part then gives
[
mu-lambda_1
=
langle b,(C-lambda_1)^{-1}bangle.
]

Set
[
t=|(C-lambda_1)^{-1}b|.
]
After choosing the phase so that (alpha>0), normalization gives
[
alpha=(1+t^2)^{-1/2}.
]
Therefore
[
|g-v|^2
=
2left(1-rac1{sqrt{1+t^2}}ight)
le t^2,
]
hence
[
|g-v|le t.
]

Finally (lambda_1lemule U). By spectral calculus for (C),
[
0<(C-lambda_1)^{-1}le(C-U)^{-1}
]
and, because both resolvents are functions of the same self-adjoint operator,
the same pointwise ordering applies to their squares. Consequently
[
|(C-lambda_1)^{-1}b|
le
|(C-U)^{-1}b|.
]
This proves the theorem.

## 4. Practical certified inverse-solve corollary

Suppose in addition that
[
C-Ugegamma I,qquad gamma>0.
]
Let (xin Qmathcal H) be any explicit approximate inverse solve and let
[
r=b-(C-U)x,qquad |r|leho.
]
Then
[
(C-U)^{-1}b=x+(C-U)^{-1}r,
]
so
[
oxed{
|g-v|
le
eta
le
|x|+rac{ho}{gamma}.
}
]
This is suitable for interval certification: (x) may be a finite
Legendre/Krylov/Chebyshev trial, (ho) can be bounded directly, and
(gamma) comes from a full-operator complement comparison.

## 5. Growing-window corollary

Let (A_L) be the localized geometric Weil operator, (v_L) the normalized
prolate candidate, and define
[
P_L=|v_Langlelangle v_L|,quad Q_L=I-P_L,quad
b_L=Q_LA_Lv_L,quad C_L=Q_LA_LQ_L.
]
Suppose along an unbounded sequence (L_n	oinfty) there are certified
numbers
[
mu_L=langle v_L,A_Lv_Langlele U_L<c_L,qquad
C_Lge c_LI,
]
and
[
oxed{
eta_L:=|(C_L-U_L)^{-1}b_L|=O(e^{-L/2}).
}
]
Then every (A_L) in the sequence has a simple isolated even ground state
(g_L), and
[
oxed{|g_L-v_L|=O(e^{-L/2}).}
]

In the normalization used in the saved growing-window transfer theorem, the
Connes–Consani–Moscovici prolate/Hermite estimate supplies
(delta_L=O(e^{-2L})). Substituting the direct distance bound above into
the existing complex-strip estimate gives, for every fixed (0<R<1/2),
[
sup_{|operatorname{Im}z|le R}|F_L(z)-Xi(z)|
=
O(e^{-(1/2-R)L})	o0.
]

At every fixed window, the operator-core result together with the
Connes–van Suijlekom real-zero theorem implies that the Fourier transform of
the simple even ground state has only real zeros. Local uniform convergence
on every closed substrip then excludes nonreal zeros of (Xi) by Hurwitz
(or local Rouché). Hence the stated growing-window hypotheses imply RH.

## 6. What now has to be proved

The global problem is reduced to a scale-dependent pair of estimates:

[
oxed{
C_Lge c_LI>U_LI
}
]
and
[
oxed{
|(C_L-U_L)^{-1}Q_LA_Lv_L|=O(e^{-L/2}).
}
]

The second estimate may be certified through the practical inverse-solve
corollary:
[
|x_L|+ho_L/gamma_L=O(e^{-L/2}),
qquad
C_L-U_Lgegamma_LI.
]

This formulation avoids requiring a positive lower bound for the ground
eigenvalue and avoids using the raw residual norm
(|Q_LA_Lv_L|), which can be dominated by high-energy components.

## 7. Current fixed-window checkpoint

At (L=1.5), the immediate numerical target remains
[
s_*=5.26	imes10^{-97},
]
with desired certified comparison inertias
[
	ext{even }(863,1,0),qquad
	ext{odd }(864,0,0).
]
Together with the certified even trial below (s_*), this establishes the
simple-even condition at the new window. The next calculation should then
construct the prolate (v_{1.5}), deflate it, and certify the
shorted-resolvent quantity (eta_{1.5}).
