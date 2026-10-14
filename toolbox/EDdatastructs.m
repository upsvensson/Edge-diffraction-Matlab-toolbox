function EDdatastructs
% EDdatastructs - Describes the data structs of the EDtoolbox.
% This function contains only documentation. Type
%   help EDdatastructs
% or just EDdatastructs, to read it. The function headers refer to this
% description instead of repeating it.
%
% OVERVIEW
% EDmain_convex takes six input structs, checks them (EDcheckinputstructs)
% and then computes a chain of internal structs. Each internal struct is
% saved in a stage file, <filestem>_<suffix>.mat, together with a hash of
% everything it depends on, so that it can be reused by later runs
% (see EDstageversions and EDrecycleresultfiles).
%
%   Struct              Made by                        Stage file
%   planedata           EDreadgeomatrices/EDreadcad,   _eddata (with edgedata)
%                       extended by EDedgeo
%   edgedata            EDedgeo                        _eddata
%   Sdata, Rdata        EDSorRgeo                      _Sdata, _Rdata
%   firstorderpathdata  EDfindconvexGApaths            _paths
%   edgetoedgedata      EDed2geo                       _ed2data
%   Hsubmatrixdata      EDinteg_submatrixstructure     _submatrixdata
%   (results)           EDmakefirstordertfs            _tf
%                       EDintegralequation_convex_tf   _qfinal, _tfinteq
%                       EDmakefirstorderirs            _ir
%                       EDfindHODpaths, EDmakeHODirs   _hodpaths, _irhod
%                       EDpropagateESIEBEM             _tfESIEBEM
%   EDres               EDmain_convex                  _EDres
%
% Each internal struct (edgedata, Sdata, Rdata, firstorderpathdata,
% edgetoedgedata, Hsubmatrixdata) has the field .EDinputdatahash, the hash
% of its stage. The next stage includes it in its own hash (EDparenthash).
%
% Sizes: nplanes, ncorners, nedges, nsources, nreceivers, nfrequencies.
%
%--------------------------------------------------------------------------
% INPUT STRUCTS
%   geoinputdata, Sinputdata, Rinputdata, envdata, controlparameters,
%   filehandlingparameters
% All fields, their default values and when they are used are defined in
% the tables in EDcheckinputstructs. After the check, every field exists;
% fields that are not relevant for the settings have the value NaN.
% EDmain_convex attaches the field .inputorigin to each of these structs in
% EDres and in the _settings file: for each field, 'user', 'default',
% 'derived' or 'notapplicable'.
%
% For piston sources, EDcheckinputstructs adds these fields to Sinputdata:
%   .coordinates            [npistons,3] the center points of the pistons
%   .pistongausscoordinates cell {npistons,1}, each [npoints,3]: the
%                           quadrature points of the piston
%   .pistongaussweights     cell {npistons,1}, each [npoints,1]
%   .ncornersperpiston      [npistons,1]
%   .pistonareas            [npistons,1] in m2
% .sourceamplitudes is expanded to [nsources,nfrequencies] (tf) or
% [nsources,1] (ir).
%
%--------------------------------------------------------------------------
% planedata
%   .corners                [ncorners,3] corner coordinates
%   .planecorners           [nplanes,nmaxcornersperplane] the corner numbers
%                           of each plane. Planes with fewer corners have
%                           zeros (or repeated corners) at the end.
%   .planeabstypes          [nplanes,6] (char) 'RIGID ', 'TOTABS' or 'SOFT  '
%   .planeeqs               [nplanes,4] plane equations [A B C D], with
%                           Ax + By + Cz = D; [A B C] is the normal vector
%   .ncornersperplanevec    [nplanes,1]
%   .minvals, .maxvals      [nplanes,3] bounding box of each plane
%   .planehasindents        [nplanes,1] number of indenting corners
%   .indentingcorners       [nplanes,nmaxcornersperplane] marks the
%                           indenting corners of each plane
%   .cornerinfrontofplane   [nplanes,ncorners] where each corner is
%                           relative to each plane (see EDinfrontofplane)
%   .modeltype              text, e.g. 'convex_ext' or 'singleplate'
% Added by EDedgeo:
%   .planeisthin            [nplanes,1] 1 if the plane is thin
%   .rearsideplane          [nplanes,1] the plane at the rear side of a thin
%                           plane (0 if not thin)
%   .canplaneobstruct       [nplanes,1] 1 if the plane can obstruct a path
%   .planeseesplane         [nplanes,nplanes]: 1 in front of (could see),
%                           0 behind, -1 aligned, -2 back-to-back
%   .reflfactors            [nplanes,1] +1 rigid, 0 totally absorbing,
%                           -1 pressure release
%
%--------------------------------------------------------------------------
% edgedata
%   .edgecorners            [nedges,2] corner numbers of each edge, in the
%                           reference order
%   .edgestartcoords        [nedges,3], .edgeendcoords [nedges,3]
%   .edgestartcoordsnudge, .edgeendcoordsnudge  [nedges,3] the same, moved
%                           slightly inwards along the edge (used in the
%                           visibility tests)
%   .edgelengthvec          [nedges,1]
%   .planesatedge           [nedges,2] the two planes of each edge; the
%                           first is the reference plane. With the right-
%                           hand thumb along the edge (start to end), the
%                           fingers come up through the reference plane.
%   .closwedangvec          [nedges,1] the closed wedge angle (radians)
%   .edgenvecs              [nedges,3] normal vector of the reference plane
%   .edgenormvecs           [nedges,3] unit vectors along the edges
%   .offedges               list of edge numbers that are not used
%   .indentingedgepairs     [n,2] connected edges that form an indent (they
%                           can never see each other)
%   .edgesatplane           [nplanes,nmaxedgesperplane] the edges of each
%                           plane
%   .edgeseesplane          [nplanes,nedges]: 1 at least one point in front
%                           of the plane, 0 behind, -1 aligned (can never
%                           see), -2 belongs to the plane
%   .edgerelatedcoordsysmatrices  [nedges,9]: for each edge, a 3x3 matrix
%                           (reshape(...,3,3)) for the edge-related
%                           coordinate system
%   .EDinputdatahash
%
%--------------------------------------------------------------------------
% Sdata, Rdata
% The input struct (Sinputdata/Rinputdata, after EDcheckinputstructs),
% with these fields added by EDSorRgeo. For Sdata the names end with s/S,
% for Rdata with r/R:
%   .visplanesfroms         [nplanes,nsources], values 0-5:
%                           0 behind a reflective or totabs plane
%                           1 aligned with a plane, but outside it
%                           2 in front of a reflective plane
%                           3 in front of a totabs plane
%                           4 inside a reflective plane
%                           5 inside a totabs plane
%   .vispartedgesfroms      [nedges,nsources] 0 to 2^nedgesubs-1: which
%                           parts of each edge are visible (all bits set =
%                           fully visible)
%   .vispartedgesfroms_start, .vispartedgesfroms_end  [nedges,nsources]
%                           the visible part of each edge, as fractions
%   .soutsidemodel          [nsources,1] 1 if (simple tests show that) the
%                           source is outside the model
%   .reftoshortlistS        [nedges,nsources] (for pistons:
%                           [nedges,nsources,npoints]) pointers into the
%                           three short lists below
%   .rSsho, .thetaSsho, .zSsho  short lists of the unique edge-related
%                           cylindrical coordinates (r, theta, z) of the
%                           sources relative to the edges
%   .pistonrotationmatrices, .reftorotationmatrices  (pistons only)
%   .EDinputdatahash
% For the ESIEBEM method, Rdata also has .nvecs and .weights for the
% surface receivers (from EDgensurfreceivers).
%
%--------------------------------------------------------------------------
% firstorderpathdata
%   .directsoundlist        [ncomp,3] [source, receiver, amplitude factor]
%                           for each visible direct sound path (the factor
%                           is 1, or 1/2 etc. for edge/corner hits)
%   .specrefllist           [nIS,3] [source, receiver, amplitude factor] for
%                           each valid first-order specular reflection
%   .specreflIScoords       [nIS,3] image source coordinates
%   .diffpaths              [nreceivers,nsources,nedges] 0/1: first-order
%                           diffraction paths
%   .edgeisactive           [nedges,1] 0/1
%   .ncomponents            [1,3] number of direct sound, specular and
%                           diffraction components
%   .EDinputdatahash
%
%--------------------------------------------------------------------------
% edgetoedgedata
%   .edgeseespartialedge    [nedges,nedges] (sparse) 0: the edges can not
%                           see each other; +n: 1 to 2^nedgesubs-1, how
%                           much they see each other, path along a plane;
%                           -n: the same, path not along a plane
%   .edgealignedwithedge    [nedges,nedges] (sparse) 0/1
%   .edgeconnectedwithedge  [nedges,nedges] (sparse) 1 if the edges share a
%                           corner (0 for aligned edge pairs)
%   .edgeconnectionangle    [nedges,nedges] (sparse) the angle at the
%                           shared corner
%   .edgeperptoplane        [nplanes,nedges] (sparse) 0/1
%   .edgeplaneperptoplane1, .edgeplaneperptoplane2  [nplanes,nedges]
%                           (sparse) 0/1
%   .reftoshortlistE        [nedges,nedges] pointers into the short lists:
%   .re1sho, .re2sho, .thetae1sho, .thetae2sho, .ze1sho, .ze2sho
%                           the edge-related cylindrical coordinates of the
%                           start and end points of one edge, relative to
%                           another edge
%   .examplecombE
%   .EDinputdatahash
%
%--------------------------------------------------------------------------
% Hsubmatrixdata  (the structure of the integral equation)
%   .edgepairlist           [nedgepairs,2] [to-edge, from-edge] for each
%                           edge pair = each block of unknowns
%   .nedgeelems             [nedges,1] number of quadrature points per edge
%   .Qvecstartnums, .Qvecendnums  [nedgepairs,1] the rows of each edge pair
%                           in the vector of unknowns (length nQ =
%                           Qvecendnums(end))
%   .Hsubmatrix             [nedgepairs,nedgepairs] the submatrix number of
%                           each block of the H-matrix (0 = no block)
%   .edgetripletlist        [nsub,3] [to-edge, via-edge, from-edge] of each
%                           submatrix
%   .listofsubmatrices      [nsub,3]
%   .submatrixcounter, .nuniquesubmatrices
%   .reftoshortlist         [nsub,1] for identical (symmetric) submatrices
%   .isthinplanetriplet, .isthinplaneedgepair
%   .quadraturematrix_pos, .quadraturematrix_weights  (sparse) the
%                           Gauss-Legendre nodes and weights; column n
%                           holds the n-point rule
%   .edgepairsharecorner, .edgepairangle
%   .edgepairsharedcornertype_from, .edgepairsharedcornertype_to
%                           0 = shared corner at the edge start, 1 = at the
%                           end, -1 = no shared corner
%   .EDinputdatahash
%
%--------------------------------------------------------------------------
% EDres  (output of EDmain_convex, and saved in <filestem>_EDres.mat)
%   .tfdirect, .tfgeom, .tfdiff, .tfinteqdiff, .tftot
%                           [nfrequencies,nreceivers] (doaddsources = 1)
%                           or [nfrequencies,nreceivers,nsources]
%   .tftot_ESIEBEM          (method 'tf_ESIEBEM')
%   .irdirect, .irgeom, .irdiff, .irhod, .irtot   (method 'ir')
%                           [nsamples,nreceivers(,nsources)]; .irhod can be
%                           a cell array with one cell per diffraction
%                           order (savealldifforders = 1)
%   .inteqdiagnostics       (method 'tf') see below
%   .timingstruct           the computation times of the stages
%   .geoinputdata ... .filehandlingparameters  the six input structs, each
%                           with .inputorigin
%   .EDversionnumber, .inputhash, .logfilename
% The irrelevant result fields are empty. With savefullEDresfile = 0, the
% file contains everything except the result matrices; EDgetresults loads
% them from the stage files.
%
% inteqdiagnostics
%   .solutionmethod, .nQ (number of unknowns), .nnzH (GMRES_fast)
%   GMRES: .gmresflag, .gmresrelres, .gmresiter, .solvetime, each
%   [nfrequencies,ncols] (ncols = nsources if the sources are separate,
%   otherwise 1); .assemblytime [nfrequencies,1] (GMRES_fast)
%   iterative: .normQ [difforder-2,nfrequencies,ncols]
%   .propdiagnostics        from EDcalcpropagatematrix
%   .qfinalfile             text: computed and saved, reused from, or
%                           not saved
%
% Peter Svensson (peter.svensson@ntnu.no) 14 Oct. 2026
%
% EDdatastructs

% 14 Oct. 2026 First version (by Claude, prompted by PS), collected from
% the headers of the functions that create the structs.

help(mfilename)
