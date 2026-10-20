# Integer exchange coordinates for native cube allocations

Status: author proof, October 4, 2026; independent review pending.
All four containment questions remain open. The elementary exchange is
established work. The result here is a spanning-tree **integer basis** for
the entire native kernel, followed by an exact statement of what this does
and does not import from algebraic statistics.

## 1. Objects and the existence question

Let B be a finite cubical subcomplex of {0,1}^n closed under faces.
Write Q(b,S) for a positioned cube with direction set S and canonical
base b, where b_i=0 for i in S; its vertices are b xor T, T subset S.
Every actual rooted cube (v,S) is a column, irrespective of whether its
coefficient in a given allocation is positive. Every positioned cube,
including a vertex, is a row. The native matrix is

\[
 M_{Q,(v,S)}=\mathbf1\{v\in Q,\ \operatorname{dir}(Q)\subseteq S\}.
\]

A unit allocation is a vector mu>=0 with M mu=1. Its integral points
have exactly one coefficient equal to one at each vertex, all others
zero. On ample classes these are precisely the ordinary representation
maps, by the established allocation dictionary. Thus Rep=Frac asks for
an integral point whenever this unit fiber is nonempty on an ample class;
it does not ask to preserve the positive support of a supplied point.

## 2. A basis made of elementary exchanges

For each actual cube Q=Q(b,S), choose a spanning tree of its vertex graph
and a root rho_Q. Orient its tree edges arbitrarily. An edge v->w has
w=v xor {i}, i in S. Define the four-column vector

\[
 g_{Q,v,w}=e_{(v,S)}-e_{(w,S)}
             -e_{(v,S\setminus\{i\})}+e_{(w,S\setminus\{i\})}.
 \tag{1}
\]

All four columns are actual: the lower two are rooted facets of Q.
This is the established native four-column exchange, not a new move.

**Theorem (integer exchange basis).** The vectors (1), one for each
chosen tree edge of each actual cube, form a basis over Z of
ker_Z M. They also form a basis of ker_R M over R. Their number is

\[
 \sum_{Q\in B}(2^{\dim Q}-1)=|I_B|-|B|.
 \tag{2}
\]

The columns (rho_Q,dir(Q)), one per cube, together with these kernel
vectors form a basis of the full coordinate lattice Z^{I_B}.

**Proof roadmap.** Process cubes from largest dimension downwards.
A tree transfers the root coefficients of one cube to its selected
root. The correction terms occur only on cubes of smaller dimension.
Tree elimination uses integers throughout; the remaining selected
columns have a unit triangular native minor.

**Proof.** First M g=0. On a row E whose directions contain i, the two
top columns both contribute or both fail to contribute, while the lower
columns both fail. If E's directions do not contain i, top and lower
columns at v have the same entry, as do those at w. These cover every
row, proving cancellation.

Let z be any integer coordinate vector. For a largest remaining cube Q,
consider all coordinates (v,dir(Q)) of its roots. Removing leaves of
the chosen tree, move the current leaf coefficient to its parent by
subtracting an integer multiple of the appropriate signed vector (1).
This clears that leaf's top coordinate. Only another top coordinate
of Q and coordinates on proper faces of Q change. Repeat until only
rho_Q remains. Process every cube in decreasing dimension. No already
processed top coordinate is reintroduced. The final residual is
supported on the selected columns (rho_Q,dir(Q)). This proves integral
spanning of the full coordinate lattice by the stated vectors.

For independence, consider a linear relation among these vectors and
the selected columns. On the top coordinates of a largest-dimensional
cube Q, it becomes a relation between the edge-incidence vectors of its
tree and e_{rho_Q}. Summing coordinates gives zero for the coefficient
of e_{rho_Q}. Leaf elimination then gives zero for every tree-edge
coefficient. Descending through all cubes proves independence. This
argument works over R as well as Z.

Finally, ordering rows and selected columns by dimension gives an
upper triangular matrix with diagonal entries one: a rooted column
contributes only to faces of its cube. It is therefore invertible over
Z. If Mz=0, applying M to the unique full-coordinate expansion of z
forces every selected-column coefficient to vanish. Exactly the
exchange vectors remain. Counting the rooted columns of each cube
gives (2). Empty complexes and zero-dimensional cubes cause no exception:
their contribution to the kernel basis is zero. QED.

No ampleness assumption is needed for this theorem.

## 3. Exact coordinates for the unit fiber

Choose rho_Q=b, the canonical base of each cube. Set all nonselected
coordinates of z^0 to zero and put

\[
 z^0_{(b,S)}=
 \sum_{T\supseteq S:\ Q(b,T)\in B,\ b\cap T=\varnothing}
       (-1)^{|T\setminus S|}.
 \tag{3}
\]

Here words are identified with their sets of one-coordinates.
Then Mz^0=1. Indeed a selected column (c,T) can cover row Q(b,S)
only if T contains S and c=b: c is zero on T, hence zero on S,
and must have the fixed outside-S coordinates b. Summing (3) over
all such T and reversing the finite sums gives, for every U above S,
the inner sum sum_{S subset T subset U}(-1)^{|U minus T|}; it is one
for U=S and zero otherwise. Face closure supplies all intermediate T.

This is the familiar signed Boolean/ordering inversion, used here as
an explicit starting point, not a new positive selector. For the ample
class C={0,1,...,6} in the three-cube, z^0 at root zero is one on the
empty support, minus one on each singleton, and one on each pair.
The negative singleton coefficients do not disappear just because
C is ample. C is ample because its seven shattered sets are exactly
the seven proper subsets of the three coordinates, all strongly
shattered. It is already represented by other allocations.

Let G have the basis vectors (1) as columns. Every real unit allocation
has a **unique** expression

\[
 \mu=z^0+G\lambda,\qquad z^0+G\lambda\ge0.
 \tag{4}
\]

Moreover mu is integral if and only if lambda is integral. The forward
implication uses the integer elimination theorem; the reverse uses the
integrality of z^0 and G. Consequently the unit existence problem is
exactly integer feasibility of the inequalities in (4), whenever their
real feasible set is nonempty. Each column of G has two +1 and two -1
entries. This sparsity is useful structure, **not** an integrality theorem.
There is no demonstrated easier ample-specific transfer at this point.

## 4. What the algebraic-statistics import supplies

Over a field k, introduce a formal variable x_j per actual rooted column.
For u>=0 write x^u=product_j x_j^{u_j}. Define the native toric ideal

\[
 I_M=\ker\bigl(k[x_j:j\in I_B]\longrightarrow k[y_Q:Q\in B],
                 \ x_j\longmapsto\prod_Q y_Q^{M_{Q,j}}\bigr).
\]

Equivalently it is generated by all binomials x^u-x^v with Mu=Mv.
One can see this without an external theorem: group the monomials of
any kernel polynomial by their equal image monomial; coefficients in
each group sum to zero, so that group is a sum of such differences.

Let J be generated only by x^{g^+}-x^{g^-} for the selected basis
exchanges, and let p=product_j x_j. Then

\[
 J:p^\infty=I_M.
 \tag{5}
\]

The left side means all polynomials f for which p^N f lies in J for
some N>=0. Thus (5) is a quadratic presentation **after saturation**,
not a claim that J already equals I_M in the polynomial ring.

**Proof.** Invert every x_j. An exchange binomial then differs by a
unit from x^g-1. Every integer kernel vector is an integer combination
of the basis vectors, and identities
x^{a+b}-1=x^a(x^b-1)+(x^a-1),
x^{-a}-1=-x^{-a}(x^a-1)
show that their Laurent ideal contains x^h-1 for every h in ker_Z M.
Hence J and I_M have identical extensions to the Laurent ring.
Contraction of J's extension is J:p^infinity. I_M is already saturated:
the image of a monomial is nonzero in the target domain, so an image
product vanishes only when its polynomial factor vanishes. Contracting
the identical extensions proves (5). QED.

For comparison, a **Markov basis** is a collection of integer kernel
moves connecting every nonempty nonnegative integer budget fiber,
with every intermediate coordinate vector nonnegative. Diaconis--Sturmfels,
Theorem 3.1, pp. 375--376, says that this requires generation of I_M as
a polynomial ideal, not merely after saturation. The hypotheses apply
to native M: columns are nonnegative integer vectors and are nonzero
because each contains its vertex row. Their Section 2.1(iv), p. 371,
also distinguishes lattice-basis moves from this stronger requirement.
Primary source:
[Algebraic Algorithms for Sampling from Conditional Distributions](https://www.math.ucdavis.edu/~deloera/MISC/LA-BIBLIO/trunk/Diaconis/Diaconis8.pdf).

Even a full Markov basis would connect existing integral points; it would
not prove that the unit fiber contains its first integral point. The
formal variables x_j in this toric comparison are generators, not the
allocation coefficients mu_j. Localization must not be interpreted as
a positivity assumption on mu or a rounding procedure for mu.

## 5. Exact controls, prior work, and the stopping point

The checker verifies cancellation and a full integer coordinate change
on Q0,Q1,Q2,Q3, the seven-word ample class, and a four-vertex path.
Kernel basis sizes are respectively 0,1,7,37,18,3. It separately checks
(3) against every positioned row.

A bounded unit-fiber diagnostic exhausts all 2,16,4096 edge orientations
of Q1,Q2,Q3. Vertex and edge equations justify this coverage of every
integer unit allocation. Direct checks of **all** native rows retain
2,12,744 points. Both the complete elementary-exchange graph and the
chosen-basis-exchange graph are connected in these three controls.
No higher-dimensional or general ample connectivity assertion follows;
these cubes already had integral witnesses. No larger census is selected.
Runs use standard-library exact arithmetic, nice-10, 256 MiB address
space and 60 CPU seconds, after checking MemAvailable above 2 GiB.

Closest established inputs are CS Round683's four-column ridge exchange
and conformal-circuit decomposition, the native pin-switch application,
and Attempt51's native marginal dictionary. The addition is the **full
integer spanning-tree basis**, its exact affine unit-fiber coordinates,
and its saturated quadratic toric presentation. None of those older
local exchange or two-endpoint decomposition statements establishes
this basis. No global priority claim is made.

**Campaign assessment.** This strengthens the algebraic description of
native allocations, but does not discharge an existence obligation.
The missing step in (4) is still a lattice point satisfying nonnegativity,
not a signed kernel generator. Retain the result as theory; do not fund
an automatic Markov-basis computation or a larger cube connectivity census
as a surrogate for integral existence. All four comparisons stay OPEN.
