# Window L=3/2: simple isolated even ground state by a monotone source shift

**James Hobday research programme · 10 November 2026**

## Status

This note records a new fixed-window theorem for the same localized geometric Weil operator on
\[
H_{3/2}=L^2([-3/2,3/2]).
\]
It changes only the auxiliary decomposition used to certify the operator. It does **not** change the Weil operator.

The result is
\[
\boxed{\lambda_1^+(3/2)<5.26\times10^{-97}
<\min\{\lambda_2^+(3/2),\lambda_1^-(3/2)\}.}
\]
Hence the full localized operator has a **simple, isolated, even ground state** at L=3/2.

No lower bound on the ground energy is claimed here.

## 1. Previously certified physical trial

The saved 192-mode even physical trial has
\[
\frac{q(f)}{\|f\|^2}
\in
[5.2550429369785162156610923235542389485061774251988910730261957707362
\times10^{-97}\ \pm 3.08\times10^{-165}],
\]
so in particular
\[
\lambda_1^+(3/2)<5.26\times10^{-97}.
\]

This trial belongs to the physical Weil form and is independent of the auxiliary source decomposition below.

## 2. Baseline source and shifted source

Write the finite-prime translation operator as P. The certified whole-space weighted-Schur estimate is
\[
\|P\|\le \frac{43}{10}.
\]

The baseline positive source was
\[
M_{c_1}=c_1I-P,\qquad c_1=\frac{24}{5}=4.8,
\]
with conservative lower floor
\[
M_{c_1}\succeq \frac12 I.
\]

Choose
\[
c_2=\frac{48059}{10000}=4.8059,\qquad
\delta=c_2-c_1=\frac{59}{10000}.
\]
At the fixed Fourier cutoff \(\Omega=768\), directed Arb arithmetic gives
\[
a(768)-c_2
>
1.2596095763767586506774792707\times10^{-5}>0.
\]
Therefore the same high-frequency omission remains nonnegative. The new source also obeys
\[
M_{c_2}\succeq \left(c_2-\frac{43}{10}\right)I
=\frac{5059}{10000}I.
\]

## 3. Source-shift shorting identity

Let E be the finite polynomial embedding and Q=I-EE*. For physical threshold s, let
\[
u=s+2^{-500},\qquad H_c=M_c-uI.
\]
Relative to E plus its orthogonal complement, write
\[
H_c=
\begin{pmatrix}
A_c & B^*\\
B & D_c
\end{pmatrix}.
\]

Because \(M_{c_2}=M_{c_1}+\delta I\),
\[
A_{c_2}=A_{c_1}+\delta I_E,\qquad
D_{c_2}=D_{c_1}+\delta I_Q,
\]
while B is unchanged.

The exact shorted finite source is
\[
S_c=A_c-B^*D_c^{-1}B.
\]
The complementary kernel term shifts by the opposite scalar amount, so \(K_c+A_c\) is independent of c. Consequently
\[
T_{\rm exact}(c_2)-T_{\rm exact}(c_1)
=
B^*(D_{c_1}^{-1}-D_{c_2}^{-1})B.
\]
Since \(D_{c_2}=D_{c_1}+\delta I\),
\[
D_{c_1}^{-1}-D_{c_2}^{-1}
=
\delta D_{c_1}^{-1}D_{c_2}^{-1}.
\]
Also
\[
D_{c_j}\preceq(c_j-u+43/10)I.
\]
Therefore
\[
\boxed{
T_{\rm exact}(c_2)-T_{\rm exact}(c_1)
\succeq
\alpha\,B^*B,
}
\]
where
\[
\alpha=
\frac{\delta}
{(c_1-u+43/10)(c_2-u+43/10)}.
\]

In the saved raw polynomial coordinates,
\[
B^*B=R:=B_{\rm src}-A_{\rm src}\Gamma^{-1}A_{\rm src}.
\]
At \(s=5.26\times10^{-97}\), Arb encloses
\[
\alpha
=
7.1201270423752550725172871253983\ldots\times10^{-5}
\]
with outward interval error below \(4\times10^{-75}\).

Thus, if \(T_0\) is the already-certified lower comparison for c1, then the PSD source-shift payment gives another rigorous lower comparison for the same Weil operator. In the accepted positive congruence this payment is
\[
\alpha J R J.
\]

## 4. Certified finite inertias

The complete calculation was rebuilt from the analytic source definitions: 3,360 integrated moments, all six \(864\times864\) interval source matrices, the low-degree source cross-check, and the \(2^{-500}\) finite-to-operator projection debit.

At
\[
s=5.26\times10^{-97},
\]
the shifted-source lower comparison has

\[
\boxed{(863+,1-,0)}\quad\text{in the even sector},
\]
and
\[
\boxed{(864+,0-,0)}\quad\text{in the odd sector}.
\]

For each parity the 800-dimensional trailing block is certified positive before the remaining 64-dimensional Schur complement is tested. The even 64-dimensional head has inertia \((63,1,0)\) with certified scaled dominance margin \(>0.2520\). The odd head is positive definite with certified scaled dominance margin \(>0.2505\).

## 5. Full-operator consequence

The inherited finite-to-operator transfer includes the full complement debit \(2^{-500}\). Hence the inertia statements imply that at most one even full-operator eigenvalue lies below s and no odd full-operator eigenvalue lies below s.

The physical trial proves at least one even eigenvalue lies below s. Counting multiplicity, exactly one full-operator eigenvalue lies below s. Therefore
\[
\boxed{
\lambda_1^+(3/2)<s<
\min\{\lambda_2^+(3/2),\lambda_1^-(3/2)\}.
}
\]

The ground eigenvalue is simple, isolated and even.

## 6. Real-zero consequence

The trigonometric core/self-adjointness argument used in the earlier fixed-window operator-core theorem is unchanged at L=3/2. Applying the Connes--van Suijlekom lowest-eigenfunction theorem to the normalized ground state \(g_{3/2}\), the entire function
\[
\Phi_{3/2}(z)=\int_{-3/2}^{3/2}g_{3/2}(u)e^{-izu}\,du
\]
has only real zeros.

This does **not** identify \(\Phi_{3/2}\) with \(\Xi\) and does not prove RH.

## 7. Reproducibility

Clean GitHub Actions run: 42791641270.

Frozen arithmetic output:
- even comparison inertia: (863,1,0);
- odd comparison inertia: (864,0,0);
- \(a(768)-c_2>1.2596095763\times10^{-5}\);
- \(\alpha>7.1201270423\times10^{-5}\);
- source matrices rebuilt from scratch and low-degree cross-check passed.

The run uses Python 3.12 and python-flint 0.8.0 with directed Arb arithmetic. No floating-point eigenvalue is accepted as a sign certificate.
