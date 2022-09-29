function EDversionhistory
% This function, EDversionhistory, prints out the EDtoolbox version
% history on the screen.
%
% Peter Svensson (peter.svensson@ntnu.no) 29 Sep. 2022

% 28 Jan 2018 First version, just after version 0.101

versionnumbers = cell(1);
versiontext = cell(1);
versiondate = cell(1);

icounter = 1;
versionnumbers{icounter} = '0.1';
versiondate{icounter} = '22 Jan 2018';
versiontext{icounter} = ['First version of the EDtoolbox. Copied from the ESIE2 toolbox ',...
'with new functions for finding first-order paths (EDfindconvexGApaths) and computing first-order tfs (EDmakefirstordertfs). ',...
'In addition, the analytic integration was implemented in EDwedge1st_fd. This version ',...
'passed EDverify tests 1,3-6 but not test 4 (the corner test).'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.101';
versiondate{icounter} = '23 Jan 2018';
versiontext{icounter} = ['A bug was fixed in EDfindconvexGApaths, which gave the wrong ',...
'amplitude of the direct sound for obscured receivers. This happened only in ',...
'special cases where there happened to be an extra edge or corner hit which was ',...
'unrelated to the obscuring plane. EDverify failed for tests 1,2,4 (did not try test 6), so ',...
'some error has sneaked into EDwedge1st_fd.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.102';
versiondate{icounter} = '24 Jan 2018';
versiontext{icounter} = ['In EDwedge1st_fd, the analytical integration was turned off for now, ',...
'because there seems to be a bug. Using the previous brute-force numerical integration ',...
'instead, EDverify passes tests 1,3-6 again.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.103';
versiondate{icounter} = '25 Jan 2018';
versiontext{icounter} = ['A new version of EDwedge1st_fd, with a much simplified analytical integration, ',...
'which seems to give accurate results. EDverify passes tests 1,3-6.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.104';
versiondate{icounter} = '25 Jan 2018';
versiontext{icounter} = ['A new version of EDsubmatrixstructure, to ensure that there are at least two points per edge.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.105';
versiondate{icounter} = '25 Jan 2018';
versiontext{icounter} = ['A new version of EDfindconvexGApaths. Previously, the direct sound was computed even if .directsound was set to 0.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.106';
versiondate{icounter} = '26 Jan 2018';
versiontext{icounter} = ['A new version of EDmain_convextf. When HOD was not calculated, an ',...
    'empty tfinteqdiff should have been saved, but that was not the case. Fixed now. Also, ',...
    'EDwedge1st_fd was fixed. The case useserialexp2 had not implemented.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.107';
versiondate{icounter} = '26 Jan 2018';
versiontext{icounter} = ['Implemented a new input parameter: Sinputdata.doallSRcombinations. ',...
                  'By setting this to zero, only S1 to R1, S2 to R2, etc will be ',...
                  'computed. Also, changed defaults for saving so that saveeddatafile ',...
                  'and saveSRdatafiles = 1, if not specified. Finally, removed the ',...
                  'addition of "results" to the outputdirectory. Now, saving happens ',...
                  'in exactly the directory you specify. If you dont specify an ',...
                  'outputdirectory, a folder called "results" is generated in the ',...
                  'input folder.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.108';
versiondate{icounter} = '2 Feb 2018';
versiontext{icounter} = ['Corrected the handling of sourceamplitudes; the frequency- ',...
                  'dependence was not handled before. Developed EDdebug to run through ',...
                  'many combinations of nsource, nreceivers, nfrequencies etc. ',...
                  'Added timingstruct as output data in the tfinteq file. Fixed ',...
                  'small bugs in EDfindconvexGApaths'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.109';
versiondate{icounter} = '9 Feb 2018';
versiontext{icounter} = ['Introduced a new parameter: controlparameters.'...
                        'skipfirstorder (default = 0).'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.110';
versiondate{icounter} = '9 Feb 2018';
versiontext{icounter} = ['Introduced inputdatahash in many functions, which makes it',...
                         'it possible to recycle result files. Uses the function DataHash ',...
                         'from Matlab Central.'];
              
icounter = icounter + 1;
versionnumbers{icounter} = '0.111';
versiondate{icounter} = '13 Feb 2018';
versiontext{icounter} = ['Introduced the function lgwt from Matlab Central, for Gauss-Legendre weights, ',...
                         'in EDdistelements. Replaced the quad_gauss,compute_gauss,r_jacobi,gauss files.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.112';
versiondate{icounter} = '14 Feb 2018';
versiontext{icounter} = ['Modified the recycling of results files so that the tf and tfinteq files always contain data, ',...
                         'instead of pointing to another file. Also changed the input struct names to geoinputdata, ',...
                         'Sinputdata, and Rinputdata. Removed the setupfile; instead stored a cell variable EDsettings ',...
                         'in the tf and tfinteq files, with all the settings stored. Introduced the possibility to suppress ',...
                         'the result recycling, with a new input parameter .suppressresultrecycling.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.113';
versiondate{icounter} = '28 Feb 2018';
versiontext{icounter} = ['Added a new main function: EDmain_convexESIEtime (and functions that it calls), '...
                         'which generates impulse responses. So far, only first-order impulse responses are generated.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.201';
versiondate{icounter} = '2 Mar 2018';
versiontext{icounter} = ['Added a new main function: EDmain_convexESIEBEM (and functions that it calls), '...
                         'which uses the ESIEBEM method. So, far only rectangular surfaces are handled.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.202';
versiondate{icounter} = '8 Mar 2018';
versiontext{icounter} = ['Modified the EDgensurfreceivers, receivers are now placed at a distance of 1e-4 rather than 1e-6.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.203';
versiondate{icounter} = '15 Mar 2018';
versiontext{icounter} = ['Modified EDedgeo (bug for cases when some plane was TOTABS or SOFT, in the cad-file) and EDSorRgeo (bug when some plane was TOTABS)'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.204';
versiondate{icounter} = '15 Mar 2018';
versiontext{icounter} = ['Modified EDedgeo and EDpoinpla (bug for cases when planes had different numbers of corners)'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.205';
versiondate{icounter} = '15 Mar 2018';
versiontext{icounter} = ['Modified EDpoinpla (bug-fix for v 0.204 introduced a new bug which was fixed not)'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.206';
versiondate{icounter} = '16 Mar 2018';
versiontext{icounter} = ['New: EDmain_convex_time. Seems correct for difforder = 2 but not higher'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.207';
versiondate{icounter} = '22 Mar 2018';
versiontext{icounter} = ['Changed variable name from hodir to irhod. Also implemented doaddsources in '...
                         'EDmakeHODirs, and savealldifforders in EDmain_convex_time.m'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.208';
versiondate{icounter} = '5 Apr 2018';
versiontext{icounter} = ['Implemented triangular-element quadrature for ESIEBEM'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.209';
versiondate{icounter} = '7 Apr 2018';
versiontext{icounter} = ['New controlparam:saveindividualfirstdiff, which means that each '...
                         'first-order diff. ir is saved individually.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.210';
versiondate{icounter} = '27 Apr 2018';
versiontext{icounter} = ['New function: EDconvertquadramodel. Also, implemented  '...
                         'general quadrilaterals for ESIEBEM.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.211';
versiondate{icounter} = '29 May 2018';
versiontext{icounter} = ['New functions: EDintegratebetaoverm and EDbetaoverm_fd. Implemented a '...
                         'first version of LCN (single sample) for the source term in EDinteg_souterm. '...
                         'Also, introduced a parameter planedata.planerefltypes so that some planes '...
                         'can be turned off for models specified in the input struct.'];
icounter = icounter + 1;
versionnumbers{icounter} = '0.211';
versiondate{icounter} = '18 June 2018';
versiontext{icounter} = ['Added geomacc as an input parameter to EDpoinpla and EDinfrontofplane'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.212';
versiondate{icounter} = '6 Febr 2019';
versiontext{icounter} = ['Modified EDfindconvexGApaths (bug: spec refl were not detected in all cases when nS>1 and nR >1'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.213';
versiondate{icounter} = '22 May 2019';
versiontext{icounter} = ['Introduced EDmain_nonconvex_time, and made adjustments to handle all-diffraction paths in non-convex geometries']; 

icounter = icounter + 1;
versionnumbers{icounter} = '0.214';
versiondate{icounter} = '3 June 2020';
versiontext{icounter} = ['Fixed a bug in many files so that folder names with spaces can be handled. Also '...
    'fixed an error in EDcheckinputstructs with sourceamplitudes.']; 

icounter = icounter + 1;
versionnumbers{icounter} = '0.215';
versiondate{icounter} = '20 Jan 2021';
versiontext{icounter} = ['Fixed a small bug: if you set some planes to non-reflecting, those planes '...
    'were still handled as reflective. Also, obstruction check for S/R-to-edge paths was turned off '...
    'for convex models.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.216';
versiondate{icounter} = '16 Mar 2021';
versiontext{icounter} = ['Fixed serious bug: EDpoinpla gave wrong hit results for '...
    'some plane geometries (corners and hit points on x-axis). Also, corner and '...
    'edge hit detection was fixed. Also fixed another bug which gave the wrong amplitude ...'...
    'for first-order diffraction when a corner was hit exactly at the corner.'...
    'Fixed a small bug in EDwedge1st_fd. Fixed EDverify cases 7 and 8.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.217';
versiondate{icounter} = '22 June 2021';
versiontext{icounter} = ['Fixed a small bug which strangely has not caused '...
    'an error earlier: EDpoinpla misbehaved when the number of points to check,'...
    'and the size of planecorners were mismatched.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.218';
versiondate{icounter} = '16 Aug 2021';
versiontext{icounter} = ['Made a small change in EDpoinpla which caused an error '...
    'for earlier Matlab versions (2018). Also fixed a bug in an error message in '...
    'EDreadcadgeomatrices.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.218';
versiondate{icounter} = '26 Aug 2021';
versiontext{icounter} = ['Made a small improvement in EDmakefirstordertfs and '...
    'EDmakefirstorderirs, and corresponding changes in EDmain_xxx. '...
    'The improvement reuses results more often. '];

icounter = icounter + 1;
versionnumbers{icounter} = '0.219';
versiondate{icounter} = '8 June 2022';
versiontext{icounter} = ['Fixed a bug in EDpoinpla that made the function '...
    'crash if planes with different numbers of corners were tested for point-in-plane. '];

icounter = icounter + 1;
versionnumbers{icounter} = '0.219';
versiondate{icounter} = '29 Sept. 2022';
versiontext{icounter} = ['Changed how the datahash was setup for EDmakeHODirs.m '...
    'The previous version created a huge datahash and failed to identify usable '...
    'existing files. '];

nversions = size(versiontext,2);
txtwidth = 80;

disp(' ')
disp('******************************************************')
disp('Version history of the EDtoolbox')
for ii = 1:nversions
   disp(' ')
   disp(['Version ',versionnumbers{ii},',   ',versiondate{ii}]) 
   txtstr = versiontext{ii};
   nchars = size(txtstr,2);
   txtstrfinished = 0;
   
   while txtstrfinished == 0
       nchars = size(txtstr,2);
       if nchars <= txtwidth
           disp(txtstr);
           txtstrfinished = 1;
       else
            iv = find(txtstr == 32);
            iv2 = find(iv<=txtwidth);
            nstop = iv(iv2(end));
            disp(txtstr(1:nstop-1));
            txtstr = txtstr(nstop+1:end);
       end
   end
end
