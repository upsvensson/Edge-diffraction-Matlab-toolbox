# EDtoolbox version history

Newest version first. The detailed history of each function is in the
comments at the top of its .m file.

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
