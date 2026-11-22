# EDtoolbox version history

Newest version first. The detailed history of each function is in the
comments at the top of its .m file.

## 0.643 (22 Nov. 2026)

New function EDslopediffractionplanes, a first building block for slope
diffraction handling. The edge source formulation is exact when an edge is
hit by the wave from a monopole, or by a diffraction wave along a rigid
surface; when an edge is hit by a diffraction wave through free space, the
missing slope diffraction makes the higher-order diffraction approximate,
and the rigid boundary condition is not fulfilled on the surfaces that
such waves hit. EDslopediffractionplanes finds those planes: a plane is
included if some part of it is visible from an edge that is active, not in
the plane (or a coplanar plane), at least partly in front of it, with a
path strictly inside the edge's open wedge that no plane obstructs. It
returns the planes, the edges that hit each of them, and sample points with
hit counts. For a convex polyhedron no plane is found; for a folded plate
(L-shape) the two inner faces are found. The function is not yet called by
EDmain_convex (no results change). New EDverify test 27.

## 0.642 (22 Nov. 2026)

Two bug fixes for non-convex models, found with solid (non-thin) L-shaped
bodies, where the total field at k = 0 should be 1 everywhere.

1. Coplanar planes (EDed2geo, version 3). When two planes are coplanar,
e.g. two faces that meet at a flat (180 degree) edge, an edge of one did
not see the edges of the other: a coplanar plane is not "in front of" a
plane, so the pair was never a candidate, and the coupling along the face
was lost. Found with a box whose top face was split into two coplanar
faces: errors of up to 21% at k = 0 (the plain box: 1%). Now the edges of
coplanar planes are treated as the edges of the same plane, and the split
box gives the same results as the plain box (within 1%).

2. A source within roundoff of a zone boundary (EDinteg_souterm, version
6). In the 'exactpole' source term, the delta-function limit took the side
of the zone boundary from the sign of sin(ny*phi/2), which for a source
angle of exactly pi (e.g. a plane wave along a face) is decided by
roundoff (+-4e-16). The side is now taken from lit3 (as for a source
exactly on the boundary) when |sin(ny*phi/2)| < 1e-12. Found with a solid
L-shaped body and a plane wave along its front and back faces: the result
toggled between correct (1%) and 40% wrong from one fold angle to the next
(175, 176, 177 degrees right; 174, 178, 179 wrong).

The convex prism, a single thin plate and the thin L-plate are unchanged.

## 0.641 (22 Nov. 2026)

Bug fix in EDfindGApaths: with specorder = 'auto', the image receivers can
reach a higher reflection order than the image sources. The list of
reflection planes of the specular paths then had fewer columns than
specorderused, and EDfindGApaths stopped with "vertical dimensions
mismatch". Found with an L-plate with a 60-degree fold and a plane wave at
70 degrees. Now the plane sequences are padded with zeros. No other
results change.

## 0.640 (22 Nov. 2026)

Two fixes for non-convex models, found with the folded plate (L-plate,
100 degrees) and a plane wave 15 degrees out of the symmetry plane, where
the total field at k = 0 scattered by up to 3% with ngauss (8-20).

1. The source term of skewed edge pairs (the angle of the to-edge varies
along it) now gets the 'accurate'/'exactpole' near-singularity treatment,
with the angle and the face factor of each to-edge node (EDinteg_souterm,
version 5). Before, such pairs were left to the legacy treatment, which
only corrects a slice with an extreme peak, so where a zone boundary of the
first edge crosses the second edge, the jump of the source term along the
second edge was smeared out over about one node spacing of the first edge.
The source term now agrees with an adaptive reference integration to about
1e-4, also next to the jump.

2. That jump along the second edge X is now integrated with nodes of X
split at the jump points, in the propagation to the receivers
(EDsegmentlayout version 2, "to-edge splits"; EDcalcpropagatematrix_block
version 5; EDintegralequation_convex_tf version 3). The source term on the
split nodes replaces the one on the ordinary nodes, as a correction that is
propagated together with the edge segments; the integral equation itself is
not changed. New inteqdiagnostics.ntosplits.

Results for the convex prism and a single thin plate are unchanged, bit for
bit. For the folded plate the scatter with ngauss went from up to 3% to
about 0.5%; the remaining part is spread over many edge pairs (probably
the H-matrix integration over the same jumps).

## 0.639 (22 Nov. 2026)

EDverify, test 25: the warning EDfindGApaths:specorder (order 2 would
contribute) is switched off during its call of EDmain_convex; the test only
uses the geometry.

Bug fix in the symmetry compression of the H matrix
(EDinteg_submatrixstructure, version 5). Edge triplets [edge3 edge2 edge1]
with the same geometry values share one submatrix. The comparison used the
angles of edge 1 and edge 3 (relative to the via edge 2) only as
|thetain - thetaout|, but beta also depends on thetain + thetaout. Found with
the folded plate: the triplet 2<-8<-2 (from the far edge of plate 1 via the
top edge of plate 2 back to the far edge, 55.6 degrees from the via edge)
shared the submatrix of 8<-2<-8 (24.4 degrees), which has the same lengths
and distances. The compensation of the path via that triplet was then only
76% correct, and two zone boundaries on the receiver circle (0.32 and 5.96
rad) had jumps of about -0.0014, independent of ngauss. Now |thetain +
thetaout - openangle| is part of the comparison too; it is the same for two
triplets that are mirror images in the wedge of the via edge, which can
still share a submatrix. Results for the convex prism and a single thin
plate are unchanged, bit for bit.

## 0.638 (22 Nov. 2026)

The split of the edge-1 integration in the propagation (0.637) now also
covers segments of edge 1 (the visible part of an edge for an image source,
with its own nodes). EDcalcpropagatematrix_block, version 4. Found with the
folded plate and 5 times as many receivers (step 0.002 rad): within about 4
degrees of the plane of plate 1, neighbouring receivers had steps of 0.001
to 0.003 (decreasing as 1/ngauss); the segments of the side edges of plate 2
were not split. With ngauss = 16, the largest second difference there went
from 0.0030 to 0.00001. Results for the convex prism and a single thin
plate are unchanged, bit for bit.

Known: two jumps of about -0.0014 remain on the circle (fixed in 0.639).

## 0.637 (22 Nov. 2026)

Two improvements of the propagation of the integral equation to the
receivers (EDcalcpropagatematrix_block, version 3), found with a
100-degree folded plate and receivers on a circle (k = 0 and k = 5), which
had isolated spikes and receiver-to-receiver noise. With ngauss = 16, the
largest second difference along the circle went from 0.044 to 0.0016 at
k = 0 (rms from 0.0043 to 0.0002), and from 0.062 to 0.0054 at k = 5,
smoother than ngauss = 48 before. The smooth deviation from 1 at k = 0
(the model residual) is unchanged. Results for the convex prism and a
single thin plate are unchanged, bit for bit; the computation time of the
circle went up by about 15%.

**Changes that affect the results.**
- The near-pole treatment of the propagation (the rows whose peak is
  narrow compared with the node spacing, with 'exactpole' or the
  Lorentzian correction) was only done for edge pairs along a face
  (acrossface_in = 1). Now also for the other edge pairs, with the four
  beta terms and the angle of each edge-1 node. Before, a receiver close to
  a zone boundary of such a pair got plain Gauss-Legendre, which gave a
  spike (for example 0.04 at a receiver 0.002 rad from a boundary).
- For an edge 1 that is skewed relative to edge 2, the receiver can be on
  a zone boundary of edge 2 for one point of edge 1; the integral over edge
  2 jumps there, which gave an erratic O(1/n) noise between neighbouring
  receivers. The edge-1 integration is now split at those points (own
  Gauss-Legendre rules, q interpolated from the edge-1 nodes), as for the
  receiver windows of edge 2. Not yet for segments of edge 1.
- New diagnostic: propinfo.propdiagnostics.nsplitreceivers.
- EDverify, test 25: specorder = 1. With the default 'auto' (since 0.635),
  its two parallel plates facing each other gave image points up to
  specordermax = 8 and a warning; the test only checks the edge-pair
  visibility.

## 0.636 (21 Nov. 2026)

Two bug fixes for non-convex models with thin plates, found with a
100-degree folded plate (two thin plates), a plane wave and k = 0, where
receivers on a circle got jumps of +0.025 and -0.014 at two zone boundaries
(the diffraction from the far edge of plate 1, reflected in plate 2 exactly
at the fold, and the plane of plate 1 beyond the fold). The jumps did not
get smaller with ngauss. The model and the incident wave are mirror
symmetric, but the solution of the integral equation was not (up to 98%
difference between mirrored edge pairs).

**Changes that affect the results (non-convex models with thin plates).**
- EDinteg_submatrixstructure (version 4): with symmetry compression, edge
  triplets with the same geometry values share a submatrix of H. Since 0.63
  (20 Nov. 2026), the thin-plane factor (thinplaneboost) depends on the
  incoming edge pair, so two triplets with the same geometry values could
  have different factors (for example 2<-6<-3 with the factor 1 got the
  submatrix of 3<-1<-2, with the factor 2). The thin-plane flag is now part
  of the comparison. For a single thin plate or an object without thin
  edges, all triplets have the same flag, so nothing changes.
- EDinteg_souterm (version 4): the legacy correction of a near-singular
  slice of the source term (one sample replaced by an accurate integral)
  used the angle of the start point of edge 3 also for skewed edges; now
  the angle (and the factor) of the node of that slice are used.

With both fixes, the solution is mirror symmetric (1e-11), the two jumps
are gone, and the rms deviation from 1 on the circle went from 0.023 to
0.019 (ngauss = 16). Results for the convex prism and a single thin plate
are unchanged, bit for bit.

## 0.635 (21 Nov. 2026)

New: `controlparameters.specorder = 'auto'`, the new default for non-convex
models (`geoinputdata.allownonconvex = 1`), and the new field
`controlparameters.specordermax` (default 8). Found with a 100-degree folded
plate (two thin plates) and a plane wave: with the earlier default
`specorder = 1`, a sector of receivers got a total field of about 0 instead
of 1. The first-order diffraction from the image source in plate 1, at the
edges of plate 2, assumes the double reflection plate 1 - plate 2, which was
not included. Diffraction from image sources of order N always assumes the
reflections of order N+1.

- `specorder = 'auto'`: the order is raised one step at a time. A new image
  source (image receiver) is kept only if its path reaches a receiver
  (source), a part of an edge, or a part of a plane (sample points). The
  tree stops growing when an order adds nothing, or at `specordermax`, with
  a warning if order `specordermax+1` would still contribute (for example
  parallel planes that face each other). Only image points that can
  contribute are kept, so 'auto' usually gives fewer image points than a
  fixed order (the folded plate: 2 image sources and 10 image receivers,
  against 6 and 24 for `specorder = 2`, with the same results).
- `specorder = N`: as before (the same paths and results, bit for bit), but
  with a warning if image points of order N+1 would see a receiver (source)
  or a part of an edge. The warning also gives the most acute concave edge,
  whose open angle theta allows up to ceil(pi/theta) reflections.
- New fields in the paths file: `.specorderused` and `.specordertruncated`.
- Convex models (`allownonconvex = 0`) are not affected.

## 0.634 (21 Nov. 2026)

Bug fix in the first-order diffraction, for convex and non-convex models:
a receiver (or a source) within about 1e-8 rad of a zone boundary, but not
exactly on it, gave NaN. Found with a continuity test where a point source
was moved across the plane of a face. Results that were finite before are
unchanged, bit for bit (checked for a convex prism).

**Changes that affect the results.**
- EDwedge1st_fd (version 5) and EDintegratebetaoverm (version 2): in the
  analytic approximation around the apex point, 1-cos(ny*fi) rounded to 0
  for |ny*fi| below about 2e-8, and the masked term gave Inf*0 = NaN even
  though the series expansion was the one used. Now the result goes
  smoothly to the value on each side of the zone boundary, and exactly on
  the boundary it is their average, as before.
- In the loop over the four terms (used when one term is exactly
  singular), a term close to 2*pi got no series expansion and was set to
  0. It now gets the same series expansion as in the vectorized branch.
- EDintegratebetaoverm: with exacthalf, all four terms were halved
  instead of term jj (the bug that was fixed in EDwedge1st_fd in 2021).

**Known limitation (found with the same test).** When a point source is
moved across the plane of a face, edges in that plane become visible or
hidden abruptly. The integral equation compensates for this, but the
quadrature converges only as about 1/n to 1/n^2 there: the total field
jumps by 1-3% at ngauss = 8 (up to 19% at 500 Hz for a V-groove), and by
less than 0.7% at ngauss = 16. Same cause as the near-grazing window in
0.633.

## 0.633 (21 Nov. 2026)

Fixes for non-convex models, found with a V-groove (a prism with a concave
dip in its top). The result for a receiver below the prism jumped between
about 1.0 and 1.44 (k = 0, where 1 is expected) when the groove depth was
changed slightly. It is now smooth (0.98-0.99). A 100-degree folded plate
also improved (1.16/0.88 to 1.03/1.01). Results for convex models are
unchanged (checked bit for bit for a convex prism).

**Changes that affect the results (non-convex models).**
- Points along an edge that lie on the reference face of another edge
  could get the angle 2*pi instead of 0, from a y coordinate of -1e-13
  (roundoff). For a non-thin wedge these are not the same direction, and
  the source term changed sign and size depending on the roundoff. Angles
  within 1e-6 of 2*pi are now set to 0 for non-thin edges
  (EDinteg_souterm, EDHkernelvalues, EDcalcpropagatematrix_block).
- An edge that starts on a face of the via-edge but continues across an
  indent of that plane now gets per-node angles (before, the angle of its
  start point for all nodes) and per-node factors: 1 on the face side, 2
  across (EDinteg_souterm, EDHkernelvalues).
- Same-plane edge pairs where one edge crosses the line of the other are
  classified as partly visible, so the node pairs are tested
  (EDedgepairvisibility).

**Known limitation.** When a reflected wave almost grazes a face (here
within about 0.15 degrees), the edge source density gets a narrow
feature that needs more Gauss points: ngauss = 12 gave errors up to
0.03, ngauss = 20 gave smooth results through the exact grazing case.

## 0.632 (21 Nov. 2026)

A distant source (a plane wave) exactly or almost along a face, e.g.
normal incidence on a box, gave NaN or very wrong results, also for convex
models. Fixed by the new source-term implementation 'exactpole' (the new
default). Results change only where a source is within about 1e-3 rad of
a zone boundary. New EDverify test 26. EDverify passes all 26 tests.

**Changes that affect the results.** New source-term implementation
'exactpole' (controlparameters.ESIEsourcetermimplementation), which is the
new default. A distant source (a plane wave) exactly along a face, e.g.
normal incidence on a box, gave NaN in the source term of the integral
equation: at the apex node, ch - 1 is below the rounding level, so the
near-singular beta term was 0/0. GMRES then failed and the integral
equation part was lost (in the test, 0 instead of 1 below the box). A
source within about 1e-6 rad of the plane of a face gave very wrong values
(11 instead of 1). 'exactpole' computes the flagged beta terms without
cancellation, integrates them with their exact complex poles (as
'exactpole' in the propagation), and uses the limit for a source exactly
on the zone boundary. This is done only for terms where the source is
within about 1e-3 rad of a zone boundary; all other source terms are the
same as with 'accurate' (bit for bit; EDverify test 1 showed that the
exact-pole treatment was less accurate for terms that are near-singular
because the to-edge point is close to the edge). 'accurate' and 'legacy'
are unchanged. EDinteg_souterm version 2.

**Input and output**

- controlparameters.ESIEsourcetermimplementation: new value 'exactpole'
  (the default). doc/EDinputfields.tex was generated again.
- New EDverify test 26: a plane wave exactly along four faces of a convex
  body.

## 0.631 (21 Nov. 2026)

Non-convex models (geoinputdata.allownonconvex = 1): partly visible
edges and edge pairs. A part of an edge that a source or receiver sees
now gets its own quadrature, instead of zero values at the nodes outside
it, so the results converge regularly with ngauss; and the visibility
between edges is classified conservatively (fully visible, hidden, or
partly visible, with node-pair visibility for the last). The results for
non-convex models change; the results for convex models are unchanged
(bit for bit in the tests). New EDverify test 25 (edge-to-edge
visibility). EDverify passes all tests.

**Partly visible edges, source side.** A source point (a source or an
image source) that sees only part [c1,c2] of an edge got a zero source
term at the edge nodes outside that part: a step inside the integration
interval, which made the convergence with ngauss erratic (O(1/n)). The
part now gets its own Gauss-Legendre rule, as a segment (new function
EDsegmentlayout). The segment's source term is taken one H-step to the
ordinary unknowns (an extra right-hand side of the integral equation)
and is propagated directly to the receivers. EDintegralequation_convex_tf
version 2.

**Partly visible edges, receiver side.** A receiver that sees only part
of an edge got zero propagation weights at the edge nodes outside that
part. The part now gets its own Gauss-Legendre rule, with q interpolated
from the edge nodes (EDcalcpropagatematrix_block, where the
receiver-dependent code of an edge pair was moved unchanged to a local
function).

Results for convex models are unchanged (bit for bit in the tests).

**Partly visible edge pairs.**

The edge-to-edge visibility test for non-convex models, which used only
the line between the two edge midpoints (so a pair was either fully
visible or hidden), is replaced by a conservative classification (new
function EDedgepairvisibility): fully visible, hidden (by a single plane),
or partly visible. For partly visible pairs, the visibility of every node
pair is computed (Hsubmatrixdata.visrowmask), and the source term and the
H-matrix rows of obstructed node pairs are set to zero. Results for convex
models are unchanged. EDed2geo version 2, EDinteg_submatrixstructure
version 3.

**Input and output**

- EDinteg_submatrixstructure has a new optional input parameter,
  planedata.
- New fields edgetoedgedata.partialvispairs, .partialvisoccluders and
  Hsubmatrixdata.visrowmask (see EDdatastructs).
- controlparameters.ESIEpropagationrefine gives an error for a model with
  partly visible edge pairs, or with edge segments.
- New EDverify test 25: the edge-to-edge visibility of three thin plates
  (EDedgepairvisibility) against the analytic answer.
- New field inteqdiagnostics.nsegments. EDinteg_souterm and
  EDcalcpropagatematrix_block accept a layout with segments
  (Hsubmatrixdata-like structs with the fields .fromnodelist and
  .segmentlimits, from EDsegmentlayout).

## 0.63 (20 Nov. 2026)

First step towards non-convex models (geoinputdata.allownonconvex = 1,
method 'tf'): specular reflections up to controlparameters.specorder,
first-order diffraction over the visible parts of the edges, and image
sources/image receivers for the source term and the propagation of the
integral equation. Edge-to-edge paths via specular reflections are not
included yet, and slope diffraction is not included, so the results for
non-convex models are approximate. Results for convex models are
unchanged, except for several connected thin plates (see below). EDverify
passes all tests.

**Changes that affect the results**

- Several connected thin plates (e.g., a folded plate): the thin-plane
  factor 2 of the integral equation now applies to every triplet whose
  incoming edge pair lies in one thin plate (before: only when all three
  edges belonged to the same plate). EDinteg_submatrixstructure version 2.
- Non-convex models (allownonconvex = 1): new stage EDfindGApaths instead
  of EDfindconvexGApaths, and image sources/receivers in the integral
  equation (EDimagepointdata).

**Input and output**

- New parameters controlparameters.specorder (default 1) and
  controlparameters.nedgesamples (default 32), used only with
  allownonconvex = 1. (specorder was earlier on the list of removed fields.)
- firstorderpathdata from EDfindGApaths has the new fields .diffpathlist,
  .Spoints, .Rpoints, .specreflplanes and .specorder (see EDdatastructs).
- allownonconvex = 1 with the 'ir' method gives an error.
- The description of geoinputdata.allownonconvex (EDcheckinputstructs)
  was updated, and doc/EDinputfields.tex was generated again.

**Bug fixes**

- EDed2geo stopped with an error for any plane with indents (variables
  without the edgedata./planedata. prefix). An unused analysis of indented
  planes is skipped. Edges of the same plane that are across an indent now
  get their computed angle (the path runs through the air).
- EDinteg_souterm: the edge points outside the visible part of an edge got
  their position set to zero (the start of the edge) instead of a zero
  source term. (Without effect so far, since the visible part was always
  the whole edge.)

## 0.62 (9 Nov. 2026)

The changes listed below were committed before the version number was
set; this entry makes them version 0.62. One change affects the results:
the new propagation method 'exactpole' is the default for
controlparameters.ESIEpropagationimplementation. It changes the
higher-order diffraction, mainly for receivers within a few degrees of
the plane of a face, where the error now converges regularly, as
1/ngauss^2. The earlier results are obtained with
ESIEpropagationimplementation = 'accurate'. Saved _tfinteq files are
computed again (the setting is part of their hash); _qfinal files are
reused.

**Changes that affect the results**

- 'exactpole' is the default for controlparameters.ESIEpropagationimplementation
  (see the description of 'exactpole' under Input and output).
  
**Input and output**

- New parameter controlparameters.nedgeelemsuser (advanced, default []):
  the number of Gauss-Legendre points on chosen edges, overriding ngauss
  and ngaussmin. Rows [maxlength npoints] set it by edge length (all edges
  not longer than maxlength); rows [corner1 corner2 npoints] set it for
  the edge between two corners (the rows of geoinputdata.corners), so the
  edge numbers need not be known. Hashed only if it is given; without it,
  the results are the same as before.
- EDinteg_submatrixstructure: the search for identical submatrices also
  compares the numbers of points on the three edges (needed when the
  numbers are given per corner pair). For all earlier input, the
  submatrices and results are unchanged.
- New parameter controlparameters.ESIEpropagationrefine (advanced,
  default 1 = not used): the solution of the integral equation is refined
  to more points per edge with the Nystrom interpolation formula
  (EDrefineQvec) before the propagation to the receivers. A factor for all
  edges, or rows [maxlength factor] (by edge length, as nedgeelemsuser).
  The solution and the _qfinal file are not changed, and the setting is
  part of the _tfinteq hash only if it is used. In convergence studies for
  a thin plate at frequency 0, most of the error came from the quadrature
  of the propagation integral, and refining by 3 on the short pieces near
  the corners and 1.5 on the middle piece gave 4-9 times smaller errors.
  timingdata (EDres.timingstruct.integralequation) has a 5th value, the
  time for the refinement, and inteqdiagnostics.nQpropagated gives the
  number of propagated unknowns when the refinement is used.

- New value 'exactpole' for controlparameters.ESIEpropagationimplementation
  (the default is 'exactpole'). For the near-singular rows of the
  propagation from the edge sources to the receivers (receivers close to
  the plane of a face), the two poles of the singular beta term are found
  exactly, by Newton iteration started from the Taylor values of
  Svensson & Calamia (2006), and the row is integrated by product
  integration with the exact pole pair as the weight function (the
  principle of singularity swap quadrature, af Klinteberg & Barnett). With
  'accurate', the Taylor model of the poles leaves a remainder that the
  Gauss-Legendre rule cannot resolve, and the error near the face planes
  converged irregularly. Cube, plane wave, 0 Hz, 629 receivers in a
  circle: with 'exactpole' the error of every receiver falls as 1/ngauss^2
  (ngauss = 16 to 160), and two-point Richardson extrapolation reaches
  a maximum error of about 1e-6 (before: about 3e-4, at the face planes).
  The extra time is about 10-20% of the propagation. Not implemented for
  thin holes, where 'accurate' is used (noted in
  propinfo.implementationnote of EDpropagateQvec). New
  propdiagnostics fields .nrows_exactpole and .nrows_fallback.

**Speed**

- EDrefineQvec (and so ESIEpropagationrefine): kernel blocks of the
  Nystrom step with identical input data (e.g. symmetric edge triplets)
  are computed only once, optionally on parallel workers, and shared. New
  output field .nkernelblocks. The results are unchanged (bitwise).
- EDHkernelvalues (used for the H-matrix and by EDrefineQvec): the
  distances l and m are computed for the node combinations of two edges
  and then picked out, and a faster test for off-edges. The results are
  unchanged (bitwise).

**Code structure (no change of the results)**

- Three new functions, with code moved from EDintegralequation_convex_tf:
  EDpropagateQvec (the propagation to the receivers), EDinteg_soutermsources
  (the source term for one or several sources) and EDHkernelvalues (the
  kernel of the H-matrix). EDrefineQvec uses the same functions instead of
  its own copies. The solution of the integral equation, which was written
  out twice, is a local function. The results are the same as before.
- EDrefineQvec: the propagation follows
  controlparameters.ESIEpropagationimplementation (before: always
  'accurate', the default). New option .parallelworkers.
- EDstageversions: the new functions are listed for
  EDintegralequation_convex_tf, so saved _qfinal and _tfinteq files are
  computed again once.

**Tests**

- EDverify test 24: the propagation rows of an edge pair across a face,
  for a receiver 0.01 rad from the plane of the face, against stored
  high-precision values (EDcalcpropagatematrix_block called directly):
  the exact-pole treatment must agree to 1e-9 (relative); the default
  must equal an explicit 'taylor'. Test 17 also changes
  ESIEpropagationimplementation to 'accurate'.
- EDverify test 23: nedgeelemsuser by length and by corner pair, and an
  error for a corner pair that is not an edge. Test 17 also changes
  nedgeelemsuser. The local function of tests 2 and 3 was renamed
  edtest2and3.

## 0.615 (20 Oct. 2026)

A bug fix for first-order diffraction with diff1method = 'GL', which gave
wrong results for all edges that are not 1 m long (the default 'accurate'
was correct). A new parameter, controlparameters.ngaussmin, and a
prototype of EDrefineQvec.

**Changes that affect the results**

- EDwedge1st_fd (bug fix): with diff1method = 'GL', the Gauss-Legendre sum
  was not multiplied by the length of the edge (the weights are for the
  interval [0,1]). The first-order diffraction was therefore only correct
  for edges of length 1: e.g. 37% wrong for a 1 x 0.5 m plate and 89% for
  a 1 x 0.95 x 0.45 m box. Present since 'GL' was introduced (24 Sep. 2026).
  Results computed with 'GL' should be recomputed. EDwedge1st_fd has
  version 4 in EDstageversions, so saved _tf files are recomputed.

**Input and output**

- New parameter controlparameters.ngaussmin (default 2): the minimum number
  of Gauss-Legendre points on every edge (odd if ngauss < 0). With the
  default, the discretization is the same as before, and saved stage files
  remain valid (ngaussmin is hashed only if it is not 2). Useful for very
  short edges, e.g. extra corners inserted near the corners of a plate,
  which otherwise get 2 points for any ngauss. EDintegralequation_convex_tf
  handles edges with more points than abs(ngauss).
- New function EDrefineQvec (a prototype): refines a computed solution
  of the integral equation (Qfinal, e.g. from a _qfinal file) to more
  Gauss-Legendre nodes on some or all edges, without solving the integral
  equation again, and optionally propagates it to the receivers. Two
  methods: 'nystrom' (the Nystrom interpolation formula, in two steps for
  the to- and from-variables) and 'lagrange' (polynomial interpolation).
  It can be used to separate the errors of the integral equation and of
  the propagation in convergence studies. EDrefineQvec_demo shows how.
  Does not change any results.

**Tests**

- EDverify test 21: 'GL' and 'accurate' first-order diffraction agree for
  edges of different lengths (fails with the old EDwedge1st_fd).
- EDverify test 22: ngaussmin gives the right number of points per edge,
  and a smaller error for a plate with short edges at the corners.
- EDverify test 17 also changes ngaussmin.

## 0.614 (18 Oct. 2026)

Two bug fixes for polygon piston sources, which change the results for
rectangular (non-square) and quadrilateral pistons, and for triangular
pistons with the default pistongaussorder = 3. Results for such pistons
from earlier versions should be recomputed. Also faster propagation in
the integral equation, and new timing lines in the log file.

**Changes that affect the results**

- Piston area (bug fix): the area of 4-sided pistons was wrong, except for
  squares; e.g. a 2 x 1 m rectangle got the area 1.66 m2 instead of 2 m2.
  The area is only used for the direct sound of a piston (in
  EDmakefirstordertfs), which was too strong by the factor
  true area/wrong area (20% for the 2 x 1 m rectangle).
- Triangle quadrature of order 3 (bug fix): the weights summed to 0.5
  instead of 1. This affected triangular pistons with pistongaussorder = 3
  (the default), where the diffraction (first order and integral equation)
  was half as strong as it should be, and surface receivers on triangular
  faces with surfacegaussorder = 3 ('tf_ESIEBEM'). The other orders were
  correct.
- EDSorRgeo: .pistongaussweights and .pistonareas are now part of the hash
  of Sdata, so the stages that use them get new hashes when they change
  (before, a _tf file computed with the wrong area would have been reused).

**Input and output**

- New function EDpolygonquadrature, which gives the quadrature points,
  weights and area of a polygon. It replaces the two copies of the
  quadrature code in EDcheckinputstructs (pistons) and EDgensurfreceivers
  (surface receivers). Apart from the two bug fixes, the points and
  weights are the same as before.
- Pistons, and faces for 'tf_ESIEBEM', with more than 4 corners now work
  for any convex polygon in any plane: the polygon is divided into
  triangles from its center point. Before, pistons had to be regular
  polygons, and the code called the missing function GEOcirclepoints.
- An error is given if the pistons have different numbers of quadrature
  points (different numbers of corners), which EDSorRgeo can not handle.

**Speed**

- New function EDcalcpropagatematrix_block: for the 'accurate'
  propagation (the default), the receivers are handled in blocks, and for
  each edge pair all the receivers of a block are computed at once
  (including the near-singular corrections). EDintegralequation_convex_tf
  uses it, with parfor over the blocks if parallelworkers > 0. The same
  results up to rounding errors. In a test with 300 receivers, the
  propagation was 30 times faster (in Octave). The 'legacy' propagation
  still uses EDcalcpropagatematrix for one receiver at a time.
- inteqdiagnostics.propdiagnostics now has totals over all the receivers
  (before: the values for the last receiver).

**Documentation**

- The edge numbering of the integral equation is described in the headers
  of EDintegralequation_convex_tf, EDinteg_souterm, EDcalcpropagatematrix,
  EDcalcpropagatematrix_block and EDinteg_submatrixstructure, and in
  EDdatastructs: edge 2 is always the diffracting edge, with
  q0(z3 <- z2) = G(z3 <- z2 <- S), q(z3 <- z2) = q0 + int H(z3 <- z2 <- z1)
  q(z2 <- z1) dz1, and p(R) = int int F(R <- z2 <- z1) q(z2 <- z1) dz1 dz2.
- New function EDdevguide: a checklist of the steps for changes to the
  toolbox (results affected or not, new or removed parameters and
  functions, new tests, new versions, the cheatsheet).

**Log file**

- New timing lines for the integral equation: the H-matrix (built once),
  the source term and the solution as averages per frequency (and per
  source, if the sources are handled separately; "N sources added"
  otherwise), and the propagation to the receivers in total. With one
  frequency, the per-frequency labels are left out. If the solution was
  taken from a _qfinal file, this is noted. timingdata (and
  EDres.timingstruct.integralequation) now has totals instead of the
  times for the first frequency.

## 0.613 (16 Oct. 2026)

A bug fix that changes the results of the integral equation (ESIE) for
geometries with edges of different lengths, and a new option for parallel
computing. Results for such geometries from earlier versions should be
recomputed.

**Changes that affect the results**

- EDcalcpropagatematrix (bug fix): the weights and node positions of an
  edge (dzvec2, ze2_re2, dzvec1) were only recomputed when the number of
  edge nodes changed. An edge with the same number of nodes as the
  previous one, but another length, got the values of the previous edge.
  This affected the propagation from the edge sources to the receivers for
  geometries with edges of different lengths that get the same number of
  nodes (ceil(ngauss*length/maxlength)); e.g. a 1 x 0.95 x 0.45 m box with
  ngauss = 8 had a 5% error in the higher-order diffraction, which did not
  decrease with ngauss. Geometries with equal edge lengths, like cubes and
  the EDverify geometries, were not affected. EDcalcpropagatematrix has
  stage version 2, so saved results are recomputed.

**Speed**

- New field `filehandlingparameters.parallelworkers` (default 0): the
  propagation to the receivers in the integral equation stage can run on
  parallel workers (Parallel Computing Toolbox), with a parfor over the
  receivers. 0 gives an ordinary serial loop; Inf uses all the workers of
  the parallel pool. The results are the same. It is not part of any hash.
  The log file says how many workers were used. Without the Parallel
  Computing Toolbox the loop runs serially, with a note on the screen and
  in the log file.

## 0.612 (15 Oct. 2026)

Speed. No changes of the computed results (only rounding errors). For a
cube with 629 receivers, 11 frequencies and 2 separate sources, the
integral equation stage went from 189 s to 17 s, and first-order
diffraction from 11.2 s to 1.7 s.

**Speed**

- EDintegralequation_convex_tf: the propagation to the receivers is done
  after the frequency loop. EDcalcpropagatematrix is called once per
  receiver (it has a new output, `pathlengths`, so that only the factor
  exp(-jk(L-Rstart)) is computed per frequency), instead of once per
  frequency, receiver and source. For a cube with 629 receivers, 11
  frequencies and 2 separate sources, the propagation time went from
  7.96 s to 0.91 s per frequency.
- EDintegralequation_convex_tf: the H-matrix is built once, before the
  frequency loop. Only the factor exp(-jkm) depends on the frequency, and
  m depends only on the column of H, so H(k) = B*diag(exp(-jk*mcol)) with
  a real, frequency-independent B. Per frequency, only the nQ exponentials
  are computed; for GMRES_fast, the sparse matrix is also built only once.
  New diagnostics field `maxcolumndistancediff` (a consistency check of
  the factorization; it should be 0).
- EDwedge1st_fd: first-order diffraction is computed for all frequencies
  at once. The adaptive integration (EDquadgk's algorithm, now with a
  vector-valued integrand, in the local function edquadgkmultifreq) runs
  once per integration interval instead of once per frequency, and beta,
  m and l are computed once per node. The GL branch also handles all
  frequencies at once. About 8 times faster for 30 frequencies; the
  results agree with the earlier ones to about 1e-14.

**Tests**

- EDverify test 18 deletes its result files from earlier EDverify runs
  first. Before, a second EDverify run in the same directory could reuse
  a complete _tfinteq file, and the test failed although the result was
  correct.

## 0.611 (15 Oct. 2026)

One change of the computed results (first-order diffraction with
`diff1method = 'GL'`), and documentation.

**Changes that affect the results**

- EDwedge1st_fd: with `diff1method = 'GL'`, receivers close to a zone
  boundary (not only exactly on it) are computed with the `'accurate'`
  method. Before, GL gave errors of up to 100% within a few degrees of a
  zone boundary. Results with `'accurate'` are unchanged. New EDverify
  test 20.

**Documentation**

- All function headers follow the same structure: the H1 line, the
  description, input and output parameters, author and date, and the call
  syntax. The local history comes after an empty line, so `help` does not
  show it. The "Uses the functions..." lists were removed;
  EDshowcalltree gives the call tree instead.
- New functions EDdatastructs (describes the internal data structs) and
  EDinputdoc (describes all input fields, with defaults). The input tables
  in EDcheckinputstructs have a description column, which EDinputdoc shows.
- EDinputdoc can write the input descriptions as LaTeX:
  `EDinputdoc('','EDinputfields.tex')`. The new cheatsheet
  (EDtoolbox_cheatsheet.pdf) is generated from this file, so it follows
  the code.
- EDverify test 19: doc/EDinputfields.tex agrees with the specification in
  EDcheckinputstructs.
- New internal function EDcheckheaders, which checks that the headers
  follow the structure.
- With the GMRES methods and difforder > 2, the log and screen text says
  that all diffraction orders are included (the same for any
  difforder > 2).

## 0.61 (14 Oct. 2026)

Large structural changes, and some changes that affect the results.

**Changes that affect the results**

- `controlparameters.ESIEsourcetermimplementation` and
  `.ESIEpropagationimplementation` have the default `'accurate'` (before:
  `'legacy'`).
- EDmakeHODirs: the source amplitudes were not applied to the higher-order
  diffraction IRs when the sources were not added (`doaddsources = 0`).
- EDwedge1st_fd: with `diff1method = 'GL'`, a receiver exactly on a zone
  boundary is computed with the `'accurate'` method instead (the GL
  quadrature is singular there).
- Result recycling: several stage files could be reused after a change that
  should have given a new calculation. This concerned changes of only the
  source amplitudes, `doaddsources`, `diff1method`/`diff1GLorder`,
  `doallSRcombinations`, the sign of `ngauss`, and (for ESIEBEM)
  `directsound`. Calculations without recycling were not affected.

**Input and output**

- EDcheckinputstructs: one table per input struct, with validators and
  defaults. Misspelt field names give an error; removed fields give a
  warning. `controlparameters.difforder` is obligatory.
- The origin of every input value (user, default, derived, not applicable)
  is stored as `.inputorigin` in each input struct, in EDres and in the
  _settings file.
- The geometry is given either as a geometry file or as matrices; a
  geometry file is read during the input check.
- The file `<filestem>_EDres.mat` is always saved, with a hash of the input
  data. `filehandlingparameters.savefullEDresfile = 0` gives a smaller file
  without the result matrices.
- EDgetresults finds the results from the file name, or from the input
  structs (via the input hash).
- The log file lists all input values with their origin, the solver
  diagnostics, and a points-per-wavelength check of the edge
  discretization. Log files are named `<filestem>_log_YYYYMMDD_HHhMMmSS.txt`,
  so reruns keep the earlier logs.

**Stage files and recycling**

- Each stage file has a hash that contains the hashes of the stages it is
  built on (chained hashes) and the version numbers, from the new function
  EDstageversions, of the stage function and the subroutines it uses. The
  overall version number is no longer part of these hashes.
- The solution of the integral equation is saved in its own stage file,
  `<filestem>_qfinal.mat`, and reused when only the receivers or the
  propagation settings are changed. New:
  `filehandlingparameters.saveinteqsolution` (default 1). Removed:
  `.saveinteqsousigs`, `.loadinteqsousigs`.
- Edge-to-edge data are no longer recomputed when only the sources or
  receivers are changed.

**Other**

- Integral equation solver: diagnostics per frequency, a pruned set of
  solution methods (GMRES_fast, GMRES_large, iterative_fixeddifforder);
  difforder 2 always uses the iterative method.
- EDverify: restructured, 18 tests, including tests of the stage hashes.
- New utility (internal): EDshowcalltree.
- load/save without eval in the stage functions.

## 0.603 (9 Oct. 2026)

Introduced singularity subtraction in EDinteg_souterm, which seems more
accurate than a previously implemented version of LCN.

## 0.602 (8 Oct. 2026)

Introduced singularity subtraction in EDcalcpropagatematrix, which is a
simpler way than ESIEBEM to handle the challenging radiation directions.

## 0.601 (25 Sep. 2026)

Introduced GL quadrature for first-order diffraction via
`controlparameters.diff1method`. Also introduced a new function EDgetGLrule
for more efficient handling of quadrature rules.

## 0.600 (18 Sep. 2026)

Introduced GMRES, in a fast version and a version for large problems. Also
introduced a top-level parameter: `controlparameters.generalmethod`.
Streamlined EDcheckinputstructs.

## 0.504 (19 April 2024)

Replaced the numerical integration in EDwedge1st_fd. Instead of the
built-in quadgk, a stripped-down version, EDquadgk, was made, where the
integrand is written directly into the function. This saves substantial
amounts of time. One test case: 64-sided disc, 10 receivers, 100
frequencies; the new calculation time was 29% of the old time.

## 0.503 (17 April 2024)

Fixed a small bug in EDmain_convex that caused an error for IR of
difforder 2 or higher, and `controlparameters.savealldifforders = 0`. Also
removed lines in EDwedge2nd that cleared variables; removing that clearing
saves lots of time for one test case: many on-axis receivers and a
polygonal circular disc (saves 56%).

## 0.502 (16 April 2024)

Saved the input data in a new _settings file, and returned those settings
in the EDres struct.

## 0.501 (25 March 2024)

Implemented the new input parameter `geoinputdata.listofedgestoskip`.

## 0.500 (29 Nov. 2023)

Changed directory structure.

## 0.402 (29 Nov. 2023)

Implemented piston sources, in the FD.

## 0.401 (30 Oct. 2023)

Cleaned up some of the input parameter lists and the EDinputdatahashes.

## 0.400 (27 Oct. 2023)

S and R structs changed, and many function calls now use entire structs as
input.

## 0.304 (14 Oct. 2023)

Small change for the free-field case: only the direct sound is saved in the
result file, and returned by the function.

## 0.303 (12 Oct. 2023)

Introduced the possibility to use a free-field case.

## 0.302 (9 Oct. 2023)

A number of smaller bug fixes.

## 0.301 (3 Oct. 2023)

Implemented the IR calculation in the single EDmain_convex. Also added the
results as an optional output parameter to EDmain_convex. Also added a
settings hash, EDsettingshash, which is saved in all results files.

## 0.300 (29 Sep. 2023)

Made a single EDmain_convex and modified many functions to "version 2",
which makes the main function much more readable. The IR calculation has
not yet been implemented in EDmain_convex.

## 0.221 (27 Sep. 2023)

Introduced the functions EDmessage and EDpostfunctext to make the main
functions more readable.

## 0.220 (28 Oct. 2022)

Fixed a bug in EDfindconvexGApaths: the previous version missed direct
sound obstacles.

## 0.219 (29 Sep. 2022, second entry with this number)

Changed how the datahash was set up for EDmakeHODirs. The previous version
created a huge datahash and failed to identify usable existing files.

## 0.219 (8 June 2022)

Fixed a bug in EDpoinpla that made the function crash if planes with
different numbers of corners were tested for point-in-plane.

## 0.218 (26 Aug. 2021, second entry with this number)

Made a small improvement in EDmakefirstordertfs and EDmakefirstorderirs,
and corresponding changes in EDmain_xxx. The improvement reuses results
more often.

## 0.218 (16 Aug. 2021)

Made a small change in EDpoinpla which caused an error for earlier Matlab
versions (2018). Also fixed a bug in an error message in
EDreadcadgeomatrices.

## 0.217 (22 June 2021)

Fixed a small bug which strangely has not caused an error earlier:
EDpoinpla misbehaved when the number of points to check and the size of
planecorners were mismatched.

## 0.216 (16 Mar. 2021)

Fixed a serious bug: EDpoinpla gave wrong hit results for some plane
geometries (corners and hit points on the x-axis). Also, corner and edge hit
detection was fixed. Also fixed another bug which gave the wrong amplitude
for first-order diffraction when a corner was hit exactly at the corner.
Fixed a small bug in EDwedge1st_fd. Fixed EDverify cases 7 and 8.

## 0.215 (20 Jan. 2021)

Fixed a small bug: if you set some planes to non-reflecting, those planes
were still handled as reflective. Also, the obstruction check for
S/R-to-edge paths was turned off for convex models.

## 0.214 (3 June 2020)

Fixed a bug in many files so that folder names with spaces can be handled.
Also fixed an error in EDcheckinputstructs with sourceamplitudes.

## 0.213 (22 May 2019)

Introduced EDmain_nonconvex_time, and made adjustments to handle
all-diffraction paths in non-convex geometries.

## 0.212 (6 Feb. 2019)

Modified EDfindconvexGApaths (bug: specular reflections were not detected
in all cases when nS > 1 and nR > 1).

## 0.211 (18 June 2018, second entry with this number)

Added geomacc as an input parameter to EDpoinpla and EDinfrontofplane.

## 0.211 (29 May 2018)

New functions: EDintegratebetaoverm and EDbetaoverm_fd. Implemented a first
version of LCN (single sample) for the source term in EDinteg_souterm. Also
introduced a parameter `planedata.planerefltypes`, so that some planes can
be turned off for models specified in the input struct.

## 0.210 (27 April 2018)

New function: EDconvertquadramodel. Also implemented general quadrilaterals
for ESIEBEM.

## 0.209 (7 April 2018)

New control parameter: `saveindividualfirstdiff`, which means that each
first-order diffraction IR is saved individually.

## 0.208 (5 April 2018)

Implemented triangular-element quadrature for ESIEBEM.

## 0.207 (22 Mar. 2018)

Changed the variable name from hodir to irhod. Also implemented
doaddsources in EDmakeHODirs, and savealldifforders in
EDmain_convex_time.

## 0.206 (16 Mar. 2018)

New: EDmain_convex_time. Seems correct for difforder = 2 but not higher.

## 0.205 (15 Mar. 2018)

Modified EDpoinpla (the bug fix for v. 0.204 introduced a new bug, which
was fixed now).

## 0.204 (15 Mar. 2018)

Modified EDedgeo and EDpoinpla (bug for cases when planes had different
numbers of corners).

## 0.203 (15 Mar. 2018)

Modified EDedgeo (bug for cases when some plane was TOTABS or SOFT in the
cad file) and EDSorRgeo (bug when some plane was TOTABS).

## 0.202 (8 Mar. 2018)

Modified EDgensurfreceivers: receivers are now placed at a distance of
1e-4 rather than 1e-6.

## 0.201 (2 Mar. 2018)

Added a new main function: EDmain_convexESIEBEM (and the functions that it
calls), which uses the ESIEBEM method. So far, only rectangular surfaces
are handled.

## 0.113 (28 Feb. 2018)

Added a new main function: EDmain_convexESIEtime (and the functions that it
calls), which generates impulse responses. So far, only first-order impulse
responses are generated.

## 0.112 (14 Feb. 2018)

Modified the recycling of results files so that the tf and tfinteq files
always contain data, instead of pointing to another file. Also changed the
input struct names to geoinputdata, Sinputdata and Rinputdata. Removed the
setup file; instead, a cell variable EDsettings is stored in the tf and
tfinteq files, with all the settings. Introduced the possibility to
suppress the result recycling, with a new input parameter
`.suppressresultrecycling`.

## 0.111 (13 Feb. 2018)

Introduced the function lgwt from Matlab Central, for Gauss-Legendre
weights, in EDdistelements. Replaced the quad_gauss, compute_gauss,
r_jacobi and gauss files.

## 0.110 (9 Feb. 2018)

Introduced inputdatahash in many functions, which makes it possible to
recycle result files. Uses the function DataHash from Matlab Central.

## 0.109 (9 Feb. 2018)

Introduced a new parameter: `controlparameters.skipfirstorder` (default 0).

## 0.108 (2 Feb. 2018)

Corrected the handling of sourceamplitudes; the frequency dependence was
not handled before. Developed EDdebug to run through many combinations of
nsources, nreceivers, nfrequencies etc. Added timingstruct as output data
in the tfinteq file. Fixed small bugs in EDfindconvexGApaths.

## 0.107 (26 Jan. 2018)

Implemented a new input parameter: `Sinputdata.doallSRcombinations`. By
setting this to zero, only S1 to R1, S2 to R2, etc. will be computed. Also
changed the defaults for saving, so that saveeddatafile and saveSRdatafiles
= 1 if not specified. Finally, removed the addition of "results" to the
outputdirectory. Now, saving happens in exactly the directory you specify.
If you don't specify an outputdirectory, a folder called "results" is
generated in the input folder.

## 0.106 (26 Jan. 2018)

A new version of EDmain_convextf. When HOD was not calculated, an empty
tfinteqdiff should have been saved, but that was not the case. Fixed now.
Also, EDwedge1st_fd was fixed: the case useserialexp2 had not been
implemented.

## 0.105 (25 Jan. 2018)

A new version of EDfindconvexGApaths. Previously, the direct sound was
computed even if `.directsound` was set to 0.

## 0.104 (25 Jan. 2018)

A new version of EDsubmatrixstructure, to ensure that there are at least
two points per edge.

## 0.103 (25 Jan. 2018)

A new version of EDwedge1st_fd, with a much simplified analytical
integration, which seems to give accurate results. EDverify passes tests
1, 3-6.

## 0.102 (24 Jan. 2018)

In EDwedge1st_fd, the analytical integration was turned off for now,
because there seems to be a bug. Using the previous brute-force numerical
integration instead, EDverify passes tests 1, 3-6 again.

## 0.101 (23 Jan. 2018)

A bug was fixed in EDfindconvexGApaths, which gave the wrong amplitude of
the direct sound for obscured receivers. This happened only in special
cases where there happened to be an extra edge or corner hit which was
unrelated to the obscuring plane. EDverify failed for tests 1, 2, 4 (test 6
was not tried), so some error has sneaked into EDwedge1st_fd.

## 0.1 (22 Jan. 2018)

First version of the EDtoolbox. Copied from the ESIE2 toolbox, with new
functions for finding first-order paths (EDfindconvexGApaths) and computing
first-order tfs (EDmakefirstordertfs). In addition, the analytic
integration was implemented in EDwedge1st_fd. This version passed EDverify
tests 1, 3-6 but not test 4 (the corner test).
