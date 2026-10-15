function passtest = EDverify(outputdirectory,runtest,showtext,plotdiagrams)
% EDverify - Runs the EDtoolbox for a number of defined tests, and compares
% the results to expected results. A logfile is written.
%
% There are two kinds of tests:
%   Tests 1-9 check the accuracy of the computed sound field, against
%   known answers or physical requirements (such as continuity).
%   Tests 10-18 check that the code paths, the files and the input
%   handling work as intended (they compare the toolbox with itself).
%   Tests 16 and 17 check the stage hashes: that the dependency lists in
%   EDstageversions agree with the code, and that a change of one input
%   parameter gives new stage files for exactly the stages that depend on
%   it.
%
% All tests use their own filestem, EDverify_test<N>, in the output
% directory, and (except tests 15 and 17, which test it) result recycling is
% switched off, so that every test really runs the present code.
%
% Input parameters:
%   outputdirectory (optional) The directory where the calculation results
%                   and the logfile will be stored. If it is not given, or
%                   empty, a window will ask for a directory.
%   runtest         (optional) A vector of 0 or 1, which for each position n
%                   tells if test n should be run. Default value: all 1.
%                   A shorter vector is padded with zeros.
%   showtext        (optional) 0 or 1; determines if the details of the
%                   results should be printed on the screen. The details
%                   are always written to the logfile. Default value: 0.
%   plotdiagrams    (optional) 0 or 1: determines if result plots will be
%                   generated (for some of the tests). Default value: 0.
%
% Output parameter:
%   passtest        A vector with -1 (fail), 0 (not run), or 1 (pass) for
%                   all the tests.
%
% Tests:
%  1. Response at cuboid surface, plane wave incidence, 0 Hz, no
%     singularity problem. HOD test.
%  2. Field continuity across zone boundaries, single edge, 100 Hz.
%     Perpendicular edge hit. Diff1 test.
%  3. Field continuity across zone boundaries, single edge, 100 Hz.
%     Skewed edge hit. Diff1 test.
%  4. Field continuity across corner zone boundaries, 100 Hz.
%  5. Field continuity for receivers close to a single edge, 0.1 Hz.
%     Diff1 test.
%  6. Replicate a non-centered internal monopole, at 0.1 Hz.
%  7. Direct sound obscuring for a corner-on hit of an octahedron.
%  8. Direct sound obscuring for an edge-on hit of a cube.
%  9. No kinks in the total field where the receivers cross the extension
%     of a cube face, but no lower-order term changes (higher-order zone
%     boundaries). Cube, source 0.1 mm from a face, 100 Hz, ngauss = 8.
% 10. The three integral equation solution methods give the same result.
% 11. Several sources computed separately (doaddsources = 0) give the same
%     results as one source at a time; for GMRES and the iterative method.
% 12. A results file without the result matrices (savefullEDresfile = 0):
%     EDgetresults gives the same EDres as EDmain_convex, found both by
%     the input data and by the filestem.
% 13. Reusing the input structs from EDres: the origin record
%     (.inputorigin) is kept for unchanged values; a changed value becomes
%     'user'.
% 14. Input checking: a misspelt field gives an error; a removed field
%     gives a warning.
% 15. Result recycling: an unchanged rerun reuses the _tfinteq file; a
%     changed ESIEpropagationimplementation gives a new one.
% 16. The dependency lists in EDstageversions agree with the code.
% 17. A change of one input parameter gives new stage files for exactly
%     the stages that depend on it.
% 18. The integral equation solution (_qfinal) is reused for new
%     receivers, and gives the same result as a run without recycling.
% 19. The input field descriptions for the cheatsheet (doc/EDinputfields.tex)
%     agree with the specification in EDcheckinputstructs.
%
% See also EDmain_convex, EDstageversions
%
% Peter Svensson (peter.svensson@ntnu.no) 15 Oct. 2026
%
% passtest = EDverify(outputdirectory,runtest,showtext,plotdiagrams);

% 15 Jan. 2018 First version
% 15 Jan. 2018 Cleaned up the EDversion/EDversionnumber confusion. Added a
% test 4, with fine-grain ZB check.
% 16 Jan. 2018 Changed the function call to EDgetversion instead of
% EDversion.
% 16 Jan. 2018 Changed EDverify into a function rather than a script.
%              Added the plotdiagrams parameter.
% 17 Jan 2018 Added test 5
% 18 Jan 2018 Added test6
% 26 Jan 2018 Removed the "results" from the loading of results file, to
% update to version 0.107.
% 31 Jan 2018 Added "results" to the default output directory.
% 5 Feb 2018  Added test 7, after that problem was discovered with v 0.108.
% Made the logfile optional, via the specification of an outputdirectory. 
% 15 Feb 2018 Added the description of test 8.
% 22 May 2019 Fixed a bug: ntests was still set to 7.
% 3 June 2020 Fixed a bug: folder names with spaces can be handled now
% 7 Oct. 2023 Adapted to the EDmain_convex of v0.300
% 22 March 2024 Removed a few text printouts which mentioned
% EDmain_convexESIE.
% 18 Apr. 2024 Reaplced some calls to GEOcalccosfi with EDcalccosfi.
% 12 Oct. 2026 Restructured (by Claude, prompted by PS): each test is a
% local function, and one loop runs the tests and writes the log. The
% results are taken from the output EDres instead of from the files.
% controlparameters.generalmethod = 'tf' instead of the old .docalctf etc.
% Result recycling is switched off (except in test 15), and each test has
% its own filestem. The logfile is written in the output directory, with
% the line ending of the computer. The outputdirectory input parameter is
% now used (it was overwritten by a directory window before). Fixed: the
% log of tests 2 and 3 used the pass value of test 4. Added tests 9-15.
% 12 Oct. 2026 After the first run: test 9 has a more robust kink measure
% (near the boundary vs. the rest of the window); the limits of tests 1
% and 10 were set from the results with the accurate implementations.
% 14 Oct. 2026 Added test 16 (the dependency lists in EDstageversions agree
% with the actual calls in the code) and test 17 (which stage files get a
% new hash when one input parameter is changed).
% 14 Oct. 2026 Test 17 also covers the _qfinal file. Added test 18 (reuse
% of the integral equation solution for new receivers).
% 14 Oct. 2026 Tests 17 and 18: the first run is done without recycling.
% Test 18: a receiver was in the plane of a cube face.
% 15 Oct. 2026 The header was updated (tests 16-18 listed).
% 15 Oct. 2026 Added test 19 (the cheatsheet input descriptions are up to
% date).

ntests = 19;

if nargin < 1 || isempty(outputdirectory)
    outputdirectory = uigetdir('', 'Select the directory for the calculation results and the logfile');
    if ~ischar(outputdirectory)
        disp('EDverify: no directory was selected.')
        passtest = zeros(1,ntests);
        return
    end
end
if nargin < 2 || isempty(runtest)
    runtest = ones(1,ntests);
end
runtest = [runtest(:).' zeros(1,ntests-numel(runtest))];
if nargin < 3
    showtext = 0;
end
if nargin < 4
    plotdiagrams = 0;
end

if exist(outputdirectory,'dir') ~= 7
    mkdir(outputdirectory)
end

passtest = zeros(1,ntests);

[EDversionnumber,changedate] = EDgetversion;
NN = num2str(EDversionnumber);
clockvec = clock;
datetimestring = sprintf('%04d%02d%02d_%02dh%02dm%02d',clockvec(1),clockvec(2),clockvec(3),...
    clockvec(4),clockvec(5),floor(clockvec(6)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The logfile, in the output directory

logfilename = fullfile(outputdirectory,['EDverify_v',NN,'_',datetimestring,'.txt']);
fid = fopen(logfilename,'w');
if fid == -1
    error('ERROR: The logfile %s could not be opened.',logfilename)
end
headerline = ['EDverify, EDtoolbox v. ',NN,' (last change on ',changedate,'), run on ',datetimestring];
fprintf(fid,'%s\n',headerline);
disp(' ')
disp(headerline)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the tests. Each test function returns a struct with the fields
%   .title  one line
%   .lines  cell array with the description and the results
%   .pass   1 or -1

testfunctions = {@edtest1,@edtest2,@edtest3,@edtest4,@edtest5,@edtest6,@edtest7,...
    @edtest8,@edtest9,@edtest10,@edtest11,@edtest12,@edtest13,@edtest14,@edtest15,...
    @edtest16,@edtest17,@edtest18,@edtest19};

for itest = 1:ntests
    if runtest(itest) ~= 1
        continue
    end
    filestem = ['EDverify_test',int2str(itest)];
    t00 = clock;
    try
        r = testfunctions{itest}(outputdirectory,filestem,plotdiagrams);
    catch ME
        r = struct('title',['(test ',int2str(itest),' stopped with an error)'],'pass',-1);
        r.lines = {['The test stopped with the error: ',ME.message]};
    end
    passtest(itest) = r.pass;
    if r.pass == 1
        verdict = 'PASSED';
    else
        verdict = 'FAILED';
    end
    resultline = sprintf('Test %d: %s ... %s (%.1f s)',itest,r.title,verdict,etime(clock,t00));

    disp(resultline)
    if showtext > 0 || r.pass ~= 1
        for ii = 1:length(r.lines)
            disp(['      ',r.lines{ii}])
        end
    end

    fprintf(fid,'\n####################################################################\n');
    fprintf(fid,'%s\n',resultline);
    for ii = 1:length(r.lines)
        fprintf(fid,'   %s\n',r.lines{ii});
    end
end

summaryline = sprintf('%d tests passed, %d failed, %d not run.',sum(passtest==1),sum(passtest==-1),sum(passtest==0));
fprintf(fid,'\n%s\n',summaryline);
fclose(fid);
disp(summaryline)
disp(['Logfile: ',logfilename])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shared helper functions

function fh = edfilehandling(outputdirectory,filestem)
% The filehandlingparameters used by all the tests: no recycling of
% existing results, no text on the screen.
fh = struct('outputdirectory',outputdirectory,'filestem',filestem,'savelogfile',1,...
    'showtext',0,'suppressresultrecycling',1);
end

function [corners,planecorners] = edcube(halfside)
% A cube, centered at the origin, with side 2*halfside.
corners = halfside*[-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1];
planecorners = [1 4 3 2; 5 6 7 8; 1 2 6 5; 3 4 8 7; 2 3 7 6; 1 5 8 4];
end

function r = edresult(title,description,passcondition,resultlines)
r = struct('title',title,'pass',-1);
if passcondition
    r.pass = 1;
end
r.lines = [description(:).' resultlines(:).'];
end

function relerr = edzbcheck(tfvalues)
% The relative difference between the value in the middle and the mean of
% the two neighbouring values (used by tests 2-4).
ivmid = ceil(length(tfvalues)/2);
meanvalue = ( abs(tfvalues(ivmid-1)) + abs(tfvalues(ivmid+1)) )/2;
relerr = abs( (meanvalue - abs(tfvalues(ivmid)))/meanvalue );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests 1-9: accuracy

function r = edtest1(outputdirectory,filestem,~)
title = 'DC response at a cuboid surface, plane wave incidence';
description = {'Expected rel. errors for receiver positions 1 (in front, at surface) and 2 (behind, at surface),',...
    'with 96 Gauss-Legendre points for the longest edge and diffraction order 20, are 3.4e-6 and 9.1e-6',...
    '(with the accurate source term and propagation, Oct. 2026; it was 1.3e-5 and 1.5e-5 before).'};

corners = [-0.2 -0.44 -0.32; 0.2 -0.44 -0.32; 0.2 0.2 -0.32; -0.2 0.2 -0.32;...
           -0.2 -0.44 0; 0.2 -0.44 0; 0.2 0.2 0; -0.2 0.2 0];
planecorners = [1 4 3 2; 5 6 7 8; 1 2 6 5; 3 4 8 7; 2 3 7 6; 1 5 8 4];
geoinputdata = struct('corners',corners,'planecorners',planecorners);
Sinputdata = struct('coordinates',1e8*[1 1 1]/sqrt(3));
Rinputdata = struct('coordinates',[0 0 0.00001; 0 0 -0.32001]);
controlparameters = struct('generalmethod','tf','frequencies',0,'difforder',20,'ngauss',96,...
    'Rstart',norm(Sinputdata.coordinates));

EDres = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,...
    edfilehandling(outputdirectory,filestem));

relerr = abs(EDres.tftot*norm(Sinputdata.coordinates) - 1);
r = edresult(title,description,relerr(1) < 5e-6 && relerr(2) < 1.4e-5,...
    {['Computed rel. errors: ',num2str(relerr(1)),' and ',num2str(relerr(2))]});
end

function r = edtest23(outputdirectory,filestem,plotdiagrams,receiverheight,hittype)
% Tests 2 and 3 differ only in the receiver height (perpendicular or
% skewed edge hit).
title = ['First-order diffraction continuity across zone boundaries, 100 Hz, ',hittype,' edge hit'];
description = {'Single edge; receivers distributed at, and very near, zone boundaries.',...
    'The value at the ZB should be very close to the mean value of the two surrounding ones (< 1e-5).'};

corners = [0 0 1; 0 0 -1; -1 0 -1; -1 0 1];
planecorners = [4 3 2 1; 1 2 3 4];
geoinputdata = struct('corners',corners,'planecorners',planecorners);
Sinputdata = struct('coordinates',[-1 1 0]/sqrt(2));
phivecd = [45.05 45.01 45.001 45.0 44.999 44.99 44.95 315.05 315.01 315.001 315.0 314.999 314.99 314.95].';
Rinputdata = struct('coordinates',[cosd(phivecd) sind(phivecd) receiverheight*ones(size(phivecd))]);
controlparameters = struct('generalmethod','tf','frequencies',100,'difforder',1,'ngauss',16,'Rstart',0);

EDres = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,...
    edfilehandling(outputdirectory,filestem));

iv1 = 8:14;
tfZBdirect = EDres.tfdirect(iv1) + EDres.tfdiff(iv1);
iv2 = 1:7;
tfZBspec = EDres.tfgeom(iv2) + EDres.tfdiff(iv2);
relerrZBdirect = edzbcheck(tfZBdirect);
relerrZBspec = edzbcheck(tfZBspec);

if plotdiagrams == 1
    figure
    plot(phivecd(iv1),abs(tfZBdirect),'-o')
    xlabel('Receiver angle   [deg.]'); ylabel('TF magnitude   [-]');
    title(['Direct sound zone boundary, ',hittype,' edge hit']);
    figure
    plot(phivecd(iv2),abs(tfZBspec),'-o')
    xlabel('Receiver angle   [deg.]'); ylabel('TF magnitude   [-]');
    title(['Specular reflection zone boundary, ',hittype,' edge hit']);
end

r = edresult(title,description,relerrZBdirect < 1e-5 && relerrZBspec < 1e-5,...
    {['Computed results at the direct sound ZB differ by: ',num2str(relerrZBdirect)],...
     ['Computed results at the specular reflection ZB differ by: ',num2str(relerrZBspec)]});
end

function r = edtest2(outputdirectory,filestem,plotdiagrams)
r = edtest23(outputdirectory,filestem,plotdiagrams,0,'perpendicular');
end

function r = edtest3(outputdirectory,filestem,plotdiagrams)
r = edtest23(outputdirectory,filestem,plotdiagrams,1,'skewed');
end

function r = edtest4(outputdirectory,filestem,~)
title = 'Continuity across corner zone boundaries, 100 Hz';
description = {'Two edges meet at 90 deg.; several receivers exactly at, and very near, zone boundaries.',...
    'The value at the ZB should be very close to the mean value of the two surrounding ones (< 1e-3).'};

corners = [0 0 -10; 0 0 1; -10 0 1; -10 0 -10];
planecorners = [4 3 2 1; 1 2 3 4];
geoinputdata = struct('corners',corners,'planecorners',planecorners,'firstcornertoskip',4);
Sinputdata = struct('coordinates',[-1 -1 0]);
Rinputdata = struct('coordinates',[0.99 -1 1.98; 0.995 -1 1.99; 1 -1 2; 1.005 -1 2.01; 1.01 -1 2.02;...
                                   0.99 1 1.98; 0.995 1 1.99; 1 1 2; 1.005 1 2.01; 1.01 1 2.02]);
controlparameters = struct('generalmethod','tf','frequencies',100,'difforder',4,'ngauss',32,'Rstart',0);

EDres = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,...
    edfilehandling(outputdirectory,filestem));

iv = 6:10;
relerrZBdirect = edzbcheck(EDres.tfdirect(iv) + EDres.tfdiff(iv) + EDres.tfinteqdiff(iv));
iv = 1:5;
relerrZBspec = edzbcheck(EDres.tfgeom(iv) + EDres.tfdiff(iv) + EDres.tfinteqdiff(iv));

r = edresult(title,description,relerrZBdirect < 1e-3 && relerrZBspec < 1e-3,...
    {['Computed results at the direct sound corner ZB differ by: ',num2str(relerrZBdirect)],...
     ['Computed results at the specular reflection corner ZB differ by: ',num2str(relerrZBspec)]});
end

function r = edtest5(outputdirectory,filestem,plotdiagrams)
title = 'First-order diffraction continuity across an edge, 0.1 Hz';
description = {'Single edge; receivers very near the edge, in front and in back.',...
    'Far source, perpendicular and skewed incidence.',...
    'Computed values in front and in back should differ by less than 2.3e-3.'};

corners = [0 0 -100; 0 0 100; -10 0 100; -10 0 -100];
planecorners = [4 3 2 1; 1 2 3 4];
geoinputdata = struct('corners',corners,'planecorners',planecorners,'firstcornertoskip',3);
soudist = 1e3;
Sinputdata = struct('coordinates',[[-1 -1 0]*soudist/sqrt(2); [-1 -1 1]*soudist/sqrt(3)]);
receivers = [-0.000001 -0.00000001 0; -0.000001 0.00000001 0];
Rinputdata = struct('coordinates',receivers);
controlparameters = struct('generalmethod','tf','frequencies',0.1,'difforder',1,'ngauss',26,'Rstart',soudist);

EDres = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,...
    edfilehandling(outputdirectory,filestem));

tftot = soudist*EDres.tftot;
tftot = [squeeze(tftot(:,:,1)) squeeze(tftot(:,:,2))];
pressurediff = [abs(tftot(:,1) - tftot(:,2)); abs(tftot(:,3) - tftot(:,4))];

if plotdiagrams == 1
    figure
    semilogx(controlparameters.frequencies,abs(tftot),'-o'); grid
    xlabel('Frequency   [Hz]'); ylabel('TF magnitude   [-]');
    title(['Sound pressure amplitude very near an edge (',num2str(norm(receivers(1,:))),' m from the edge)']);
end

r = edresult(title,description,max(pressurediff) < 2.3e-3,...
    {['Computed results differ by max: ',num2str(max(pressurediff))]});
end

function r = edtest6(outputdirectory,filestem,plotdiagrams)
title = 'Replicate an internal, non-centered monopole, 0.1 Hz';
description = {'Cube with an internal, non-centered monopole; the entire cube surface gets sources.',...
    'The radiated field should be within [0.996,1.008] around a circle of receivers.'};

[corners,planecorners] = edcube(0.5);
geoinputdata = struct('corners',corners,'planecorners',planecorners);
planedata = EDreadgeomatrices(corners,planecorners);

controlparameters = struct('generalmethod','tf','frequencies',0.1,'ngauss',32,'difforder',15);
envdata = struct('cair',344,'rhoair',1.21);

% Sources: patches on all faces, with the normal velocity of the monopole
patchpoints = []; patchcorners = []; patchareas = []; patchnvecs = [];
n1 = 6; n2 = 6;
planenvecs = planedata.planeeqs(:,1:3);
for ii = 1:6
    [patchcornercoords,patchmidcoords,areas] = EDdivrect(corners(planecorners(ii,1),:),corners(planecorners(ii,2),:),...
        corners(planecorners(ii,3),:),corners(planecorners(ii,4),:),n1,n2);
    nvec = planenvecs(ii,:);
    patchmidcoords = patchmidcoords + 0.0001*nvec(ones(n1*n2,1),:);
    patchpoints = [patchpoints; patchmidcoords]; %#ok<AGROW>
    patchcorners = [patchcorners; patchcornercoords]; %#ok<AGROW>
    patchareas = [patchareas; areas]; %#ok<AGROW>
    patchnvecs = [patchnvecs; nvec(ones(n1*n2,1),:)]; %#ok<AGROW>
end

k = 2*pi*controlparameters.frequencies/envdata.cair;
internalmonopole = [0.3 0.3 0.3];
quadweights = [64 25 25 25 25 40 40 40 40]/81/4;
quadshift = sqrt(3)/sqrt(5);

% Gauss point 1 = center point
[distances,cosfi] = EDcalccosfi(internalmonopole,patchpoints,patchnvecs);
pmono_patch = exp(-1i*k*distances)./distances;
umono_radial_patch = pmono_patch/envdata.cair/envdata.rhoair.*(1 + 1./(1i*k*distances));
umono_normal_patch = quadweights(1)*umono_radial_patch.*cosfi;
% Gauss points 2-5 = from center point, towards corners
for ii = 1:4
    patchquadraturepoints = patchpoints + quadshift*(patchcorners(:,(1:3) + 3*(ii-1)) - patchpoints);
    [distances,cosfi] = EDcalccosfi(internalmonopole,patchquadraturepoints,patchnvecs);
    pmono_patch = exp(-1i*k*distances)./distances;
    umono_radial_patch = pmono_patch/envdata.cair/envdata.rhoair.*(1 + 1./(1i*k*distances));
    umono_normal_patch = umono_normal_patch + quadweights(ii+1)*umono_radial_patch.*cosfi;
end
% Gauss points 6-9 = from center point, towards edge midpoints
for ii = 1:4
    if ii < 4
        patchedgemidpoint = 0.5*( patchcorners(:,(1:3) + 3*(ii-1)) + patchcorners(:,(1:3) + 3*ii) );
    else
        patchedgemidpoint = 0.5*( patchcorners(:,10:12) + patchcorners(:,1:3) );
    end
    patchquadraturepoints = patchpoints + quadshift*(patchedgemidpoint - patchpoints);
    [distances,cosfi] = EDcalccosfi(internalmonopole,patchquadraturepoints,patchnvecs);
    pmono_patch = exp(-1i*k*distances)./distances;
    umono_radial_patch = pmono_patch/envdata.cair/envdata.rhoair.*(1 + 1./(1i*k*distances));
    umono_normal_patch = umono_normal_patch + quadweights(ii+5)*umono_radial_patch.*cosfi;
end

Sinputdata = struct('coordinates',patchpoints,'doaddsources',1);
Sinputdata.sourceamplitudes = patchareas.*umono_normal_patch.*1i*2*pi*controlparameters.frequencies*envdata.rhoair/4/pi;

rdist = 3;
phivec = (0:0.1:2*pi).';
receivers = rdist*[cos(phivec) sin(phivec) 0*phivec];
receivers = receivers + internalmonopole(ones(size(receivers,1),1),:);
Rinputdata = struct('coordinates',receivers);

EDres = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,...
    edfilehandling(outputdirectory,filestem));

maxval = rdist*max(abs(EDres.tftot));
minval = rdist*min(abs(EDres.tftot));

if plotdiagrams == 1
    figure
    plot(phivec*180/pi,rdist*abs(EDres.tftot.'),'-o'); grid
    xlabel('Receiver angle   [deg.]'); ylabel('Sound pressure amp., re. free-field   [-]');
end

r = edresult(title,description,maxval < 1.008 && minval > 0.996,...
    {['Computed results are within [',num2str(minval),',',num2str(maxval),']']});
end

function r = edtest78(outputdirectory,filestem,corners,planecorners,source,receiver,title)
% Tests 7 and 8: the direct sound must be exactly zero.
description = {'The direct sound should be zero.'};
geoinputdata = struct('corners',corners,'planecorners',planecorners);
Sinputdata = struct('coordinates',source);
Rinputdata = struct('coordinates',receiver);
controlparameters = struct('generalmethod','tf','frequencies',100,'difforder',0,'ngauss',16,'Rstart',0);

EDres = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,...
    edfilehandling(outputdirectory,filestem));

r = edresult(title,description,abs(EDres.tfdirect) == 0,...
    {['The computed direct sound had the value ',num2str(abs(EDres.tfdirect))]});
end

function r = edtest7(outputdirectory,filestem,~)
c = 0.707106781;
corners = [c 0 0; 0 c 0; -c 0 0; 0 -c 0; 0 0 c; 0 0 -c];
planecorners = [1 2 5; 2 3 5; 3 4 5; 4 1 5; 2 1 6; 3 2 6; 4 3 6; 1 4 6];
r = edtest78(outputdirectory,filestem,corners,planecorners,[2 0 0],[-2 0 0],...
    'Direct sound obscured for a corner-on hit, octahedron');
end

function r = edtest8(outputdirectory,filestem,~)
[corners,planecorners] = edcube(0.5);
r = edtest78(outputdirectory,filestem,corners,planecorners,[1 1 0],[-1 -1 0],...
    'Direct sound obscured for an edge-on hit, cube');
end

function r = edtest9(outputdirectory,filestem,plotdiagrams)
title = 'No kinks at higher-order zone boundaries, cube, 100 Hz';
maxrelstep = 1e-3;
description = {'1 m cube, source 0.1 mm from the center of the top face, receivers on an arc of radius 2 m,',...
    'in two windows around 14.48 deg. and 104.48 deg., where the receivers cross the extension of a',...
    'side face and the bottom face. The true total field is smooth there (it is smooth everywhere',...
    'outside the object), so a step at the boundary is a numerical error. Measure of the step: the',...
    'largest second difference of abs(tftot) within 1 deg. of the boundary, minus the largest one in',...
    'the rest of the window, relative to max abs(tftot) in the window. With the legacy propagation',...
    'it was estimated to about 1.6e-2; with the accurate one it was 2.9e-4 and 3.5e-4 (Oct. 2026).',...
    ['It should be below ',num2str(maxrelstep),'.'],...
    'Also shown: the ratio of the two second differences (about 1 for a perfectly smooth curve).'};

[corners,planecorners] = edcube(0.5);
geoinputdata = struct('corners',corners,'planecorners',planecorners);
Sinputdata = struct('coordinates',[0 0 0.5001]);
thetavec1 = (10:0.25:20).';
thetavec2 = (100:0.25:110).';
thetavec = [thetavec1; thetavec2];
Rinputdata = struct('coordinates',2*[sind(thetavec) 0*thetavec cosd(thetavec)]);
controlparameters = struct('generalmethod','tf','frequencies',100,'ngauss',8,'difforder',10);

EDres = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,...
    edfilehandling(outputdirectory,filestem));

n1 = length(thetavec1);
ratios = zeros(1,2);
relsteps = zeros(1,2);
windows = {1:n1, n1+1:length(thetavec)};
boundaryangles = [asind(0.25) acosd(-0.25)];     % 14.48 and 104.48 deg.
for ii = 1:2
    thetawindow = thetavec(windows{ii});
    d2 = abs(diff(abs(EDres.tftot(windows{ii})),2));
    thetad2 = thetawindow(2:end-1);              % the angle of each second difference
    nearboundary = abs(thetad2 - boundaryangles(ii)) < 1;
    ratios(ii) = max(d2(nearboundary))/max(d2(~nearboundary));
    relsteps(ii) = max(0,max(d2(nearboundary)) - max(d2(~nearboundary)))/max(abs(EDres.tftot(windows{ii})));
end

if plotdiagrams == 1
    figure
    plot(thetavec1,abs(EDres.tftot(windows{1})),'-o',thetavec2-90,abs(EDres.tftot(windows{2})),'-*'); grid
    xlabel('Receiver angle (window 2 shifted by -90 deg.)   [deg.]'); ylabel('abs(tftot)   [-]');
end

r = edresult(title,description,all(relsteps < maxrelstep),...
    {['Relative step around 14.48 deg.: ',num2str(relsteps(1),3),', around 104.48 deg.: ',num2str(relsteps(2),3)],...
     ['Ratio of second differences around 14.48 deg.: ',num2str(ratios(1),3),', around 104.48 deg.: ',num2str(ratios(2),3)]});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests 10-15: code paths, files and input handling

function [geoinputdata,Sinputdata,Rinputdata,controlparameters] = edsmallcube
% A small, quick case used by tests 10-15.
[corners,planecorners] = edcube(0.5);
geoinputdata = struct('corners',corners,'planecorners',planecorners);
Sinputdata = struct('coordinates',[1 0.7 0.3]);
Rinputdata = struct('coordinates',[-1.5 -1 0; 2 0.2 -0.4; 0.1 1.8 1.2]);
controlparameters = struct('generalmethod','tf','frequencies',[100 200],'ngauss',8,'difforder',10);
end

function r = edtest10(outputdirectory,filestem,~)
title = 'The three integral equation solution methods agree';
tolgmres = 1e-9;
toliterative = 1e-8;
description = {'1 m cube, 2 frequencies, 3 receivers. GMRES_fast and GMRES_large should agree to',...
    [num2str(tolgmres),' (relative to max abs(tftot)); iterative_fixeddifforder with difforder 30 should agree'],...
    ['with GMRES to ',num2str(toliterative),'. (Oct. 2026: 6e-17 and 1e-11.)']};

[geoinputdata,Sinputdata,Rinputdata,controlparameters] = edsmallcube;
methods = {'GMRES_fast','GMRES_large','iterative_fixeddifforder'};
tftot = cell(1,3);
for ii = 1:3
    cp = controlparameters;
    cp.inteqsolutionmethod = methods{ii};
    if ii == 3
        cp.difforder = 30;
    end
    EDres = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,cp,...
        edfilehandling(outputdirectory,[filestem,'_',methods{ii}]));
    tftot{ii} = EDres.tftot;
end
scale = max(abs(tftot{1}(:)));
diffgmres = max(abs(tftot{2}(:)-tftot{1}(:)))/scale;
diffiterative = max(abs(tftot{3}(:)-tftot{1}(:)))/scale;

r = edresult(title,description,diffgmres < tolgmres && diffiterative < toliterative,...
    {['GMRES_large vs GMRES_fast: ',num2str(diffgmres,3)],['iterative vs GMRES_fast: ',num2str(diffiterative,3)]});
end

function r = edtest11(outputdirectory,filestem,~)
title = 'Several sources computed separately = one source at a time';
tol = 1e-9;
description = {'1 m cube, 2 sources, doaddsources = 0, compared with two runs with one source each,',...
    ['for GMRES_fast and iterative_fixeddifforder. The results should agree to ',num2str(tol),'.']};

[geoinputdata,~,Rinputdata,controlparameters] = edsmallcube;
sources = [1 0.7 0.3; -0.4 -1.2 0.8];
methods = {'GMRES_fast','iterative_fixeddifforder'};
maxdiff = zeros(1,2);
for imethod = 1:2
    cp = controlparameters;
    cp.inteqsolutionmethod = methods{imethod};
    Sinputdata = struct('coordinates',sources,'doaddsources',0);
    EDresboth = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,cp,...
        edfilehandling(outputdirectory,[filestem,'_',methods{imethod},'_both']));
    for isou = 1:2
        EDresone = EDmain_convex(geoinputdata,struct('coordinates',sources(isou,:)),Rinputdata,struct,cp,...
            edfilehandling(outputdirectory,[filestem,'_',methods{imethod},'_',int2str(isou)]));
        scale = max(abs(EDresone.tftot(:)));
        onefromboth = EDresboth.tftot(:,:,isou);
        maxdiff(imethod) = max(maxdiff(imethod),max(abs(onefromboth(:)-EDresone.tftot(:)))/scale);
    end
end

r = edresult(title,description,all(maxdiff < tol),...
    {['Largest relative difference, GMRES_fast: ',num2str(maxdiff(1),3),', iterative: ',num2str(maxdiff(2),3)]});
end

function r = edtest12(outputdirectory,filestem,~)
title = 'Results file without result matrices, and EDgetresults';
description = {'1 m cube, savefullEDresfile = 0. EDgetresults, called with the input structs and with',...
    'the filestem, should both give the same EDres as EDmain_convex, and the _EDres file should',...
    'not contain the result matrices.'};

[geoinputdata,Sinputdata,Rinputdata,controlparameters] = edsmallcube;
fh = edfilehandling(outputdirectory,filestem);
fh.savefullEDresfile = 0;

EDresmain = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,fh);
EDresbyinput = EDgetresults(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,fh);
EDresbyname = EDgetresults(outputdirectory,filestem);

filecontents = load(fullfile(outputdirectory,[filestem,'_EDres.mat']),'EDres');
iscompact = ~isfield(filecontents.EDres,'tftot');
samebyinput = isequaln(EDresmain,EDresbyinput);
samebyname = isequaln(EDresmain,EDresbyname);

r = edresult(title,description,iscompact && samebyinput && samebyname,...
    {['File without result matrices: ',mat2str(iscompact)],...
     ['Same EDres via the input structs: ',mat2str(samebyinput),', via the filestem: ',mat2str(samebyname)]});
end

function r = edtest13(outputdirectory,filestem,~)
title = 'Reused input structs keep their origin record';
description = {'1 m cube. Run A with minimal input; run B with the six structs from EDres of run A;',...
    'run C as B, but with ngauss changed. The origin records of A and B should be identical,',...
    'and in C, ngauss should be ''user'' while discretizationtype stays ''default''.'};

[geoinputdata,Sinputdata,Rinputdata,controlparameters] = edsmallcube;
controlparameters = rmfield(controlparameters,'ngauss');   % so that ngauss gets its default value

EDresA = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,...
    edfilehandling(outputdirectory,[filestem,'_A']));
fhB = EDresA.filehandlingparameters;
fhB.filestem = [filestem,'_B'];
EDresB = EDmain_convex(EDresA.geoinputdata,EDresA.Sinputdata,EDresA.Rinputdata,EDresA.envdata,...
    EDresA.controlparameters,fhB);
cpC = EDresA.controlparameters;
cpC.ngauss = 10;
fhC = EDresA.filehandlingparameters;
fhC.filestem = [filestem,'_C'];
EDresC = EDmain_convex(EDresA.geoinputdata,EDresA.Sinputdata,EDresA.Rinputdata,EDresA.envdata,cpC,fhC);

structnames = {'geoinputdata','Sinputdata','Rinputdata','envdata','controlparameters'};
sameorigin = true;
for ii = 1:length(structnames)
    sameorigin = sameorigin && isequal(EDresA.(structnames{ii}).inputorigin,EDresB.(structnames{ii}).inputorigin);
end
sameresult = isequaln(EDresA.tftot,EDresB.tftot);
originsC = EDresC.controlparameters.inputorigin;
correctC = strcmp(EDresA.controlparameters.inputorigin.ngauss,'default') && strcmp(originsC.ngauss,'user') && ...
    strcmp(originsC.discretizationtype,'default');

r = edresult(title,description,sameorigin && sameresult && correctC,...
    {['Origin records A and B identical: ',mat2str(sameorigin),', results identical: ',mat2str(sameresult)],...
     ['Run C: ngauss is ''',originsC.ngauss,''', discretizationtype is ''',originsC.discretizationtype,'''']});
end

function r = edtest14(outputdirectory,filestem,~)
title = 'Input checking: misspelt and removed fields';
description = {'A misspelt field (controlparameters.diforder) should stop with the error',...
    'EDinputcheck:unknownparam. A removed field (filehandlingparameters.savecadgeofile) should',...
    'give the warning EDinputcheck:removedparam, and the calculation should run.'};

[geoinputdata,Sinputdata,Rinputdata,controlparameters] = edsmallcube;
controlparameters.difforder = 0;

cpwrong = controlparameters;
cpwrong.diforder = 5;
errorid = '';
try
    EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,cpwrong,...
        edfilehandling(outputdirectory,[filestem,'_misspelt']));
catch ME
    errorid = ME.identifier;
end

fh = edfilehandling(outputdirectory,[filestem,'_removed']);
fh.savecadgeofile = 0;
lastwarn('','');
warningstate = warning('off','EDinputcheck:removedparam');
EDres = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,fh);
warning(warningstate);
[~,warningid] = lastwarn;

r = edresult(title,description,strcmp(errorid,'EDinputcheck:unknownparam') && ...
    strcmp(warningid,'EDinputcheck:removedparam') && ~isempty(EDres.tftot),...
    {['Error identifier for the misspelt field: ''',errorid,''''],...
     ['Warning identifier for the removed field: ''',warningid,'''']});
end

function r = edtest15(outputdirectory,filestem,~)
title = 'Result recycling of the _tfinteq file';
description = {'1 m cube, with result recycling switched on. An unchanged rerun should reuse the',...
    '_tfinteq file (it is not rewritten), and a changed ESIEpropagationimplementation should give a',...
    'new file (with another hash).'};

[geoinputdata,Sinputdata,~,controlparameters] = edsmallcube;
Rinputdata = struct('coordinates',[-1.7 0.3 0.2; 0.4 -2.1 0.9]);   % not used by the other tests
fh = edfilehandling(outputdirectory,filestem);
fh.suppressresultrecycling = 0;
tfinteqfile = fullfile(outputdirectory,[filestem,'_tfinteq.mat']);
if exist(tfinteqfile,'file') == 2
    delete(tfinteqfile)
end

controlparameters.ESIEpropagationimplementation = 'accurate';
EDres1 = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,fh);
fileinfo1 = dir(tfinteqfile);
hash1 = load(tfinteqfile,'EDinputdatahash');

pause(1.5)   % so that a rewritten file would get another time stamp
EDres2 = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,fh);
fileinfo2 = dir(tfinteqfile);
reused = (fileinfo2.datenum == fileinfo1.datenum) && isequaln(EDres1.tftot,EDres2.tftot);

controlparameters.ESIEpropagationimplementation = 'legacy';
EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,fh);
hash3 = load(tfinteqfile,'EDinputdatahash');
recomputed = ~strcmp(hash1.EDinputdatahash,hash3.EDinputdatahash);

r = edresult(title,description,reused && recomputed,...
    {['Unchanged rerun reused the file: ',mat2str(reused)],...
     ['Changed propagation implementation gave a new file: ',mat2str(recomputed)]});
end

function r = edtest16(~,~,~)
title = 'The dependency lists in EDstageversions agree with the code';
description = {'For each stage function in EDstageversions, all the functions it uses (directly or',...
    'indirectly, found by matlab.codetools.requiredFilesAndProducts) must be either in its list of',...
    'subroutines or in the ignored list. Listed functions that are not called are reported, but',...
    'are not an error (e.g., a mex file that does not exist on this computer).'};

info = EDstageversions;
stagenames = fieldnames(info.stagesubroutines);
allok = true;
resultlines = {};
for ii = 1:numel(stagenames)
    stagename = stagenames{ii};
    try
        EDstageversions(stagename);    % checks that all listed functions have a version
    catch ME
        allok = false;
        resultlines{end+1} = [stagename,': ',ME.message]; %#ok<AGROW>
        continue
    end
    stagefile = which(stagename);
    if isempty(stagefile)
        allok = false;
        resultlines{end+1} = [stagename,': the function can not be found on the Matlab path']; %#ok<AGROW>
        continue
    end
    filelist = matlab.codetools.requiredFilesAndProducts(stagefile);
    actual = {};
    for jj = 1:numel(filelist)
        if strcmp(filelist{jj},stagefile) || strncmp(filelist{jj},matlabroot,length(matlabroot))
            continue
        end
        [~,name] = fileparts(filelist{jj});
        actual{end+1} = name; %#ok<AGROW>
    end
    actual = unique(actual);
    listed = info.stagesubroutines.(stagename);
    missing = setdiff(actual,[listed(:).',info.ignored(:).']);
    extra = setdiff(listed,actual);
    if isempty(missing)
        resultlines{end+1} = [stagename,': OK (',int2str(numel(listed)),' subroutines)']; %#ok<AGROW>
    else
        allok = false;
        resultlines{end+1} = [stagename,': MISSING in EDstageversions: ',strjoin(missing,', ')]; %#ok<AGROW>
    end
    if ~isempty(extra)
        resultlines{end+1} = ['   (listed, but not called: ',strjoin(extra,', '),')']; %#ok<AGROW>
    end
end

r = edresult(title,description,allok,resultlines);
end

function r = edtest19(~,~,~)
title = 'The cheatsheet input descriptions are up to date';
description = {'doc/EDinputfields.tex (input by doc/EDtoolbox_cheatsheet.tex) must be the same as a',...
    'new file from EDinputdoc, apart from the comment lines. If not, run',...
    '   EDinputdoc('''',''<EDtoolbox>/doc/EDinputfields.tex'')',...
    'and compile EDtoolbox_cheatsheet.tex again (the PDF is in the top folder).'};

% The doc folder is next to the toolbox folder (which has EDinputdoc.m)
toolboxfolder = fileparts(which('EDinputdoc'));
committedfile = fullfile(toolboxfolder,'..','doc','EDinputfields.tex');
if exist(committedfile,'file') ~= 2
    % Not an error in the toolbox: the doc folder is simply not there
    resultlines = {['Not checked: the file ',committedfile,' was not found.']};
    r = edresult(title,description,true,resultlines);
    return
end

newfile = [tempname,'.tex'];
cleanupobj = onCleanup(@() delete(newfile));
EDinputdoc('',newfile);
committedlines = edtexlines(committedfile);
newlines = edtexlines(newfile);

if isequal(committedlines,newlines)
    resultlines = {sprintf('OK: %d lines agree.',numel(newlines))};
    r = edresult(title,description,true,resultlines);
    return
end

resultlines = {'doc/EDinputfields.tex is out of date. The first differences:'};
nshown = 0;
for ii = 1:max(numel(committedlines),numel(newlines))
    if ii > numel(committedlines) || ii > numel(newlines) || ~strcmp(committedlines{ii},newlines{ii})
        if ii <= numel(committedlines)
            resultlines{end+1} = ['   file: ',committedlines{ii}]; %#ok<AGROW>
        end
        if ii <= numel(newlines)
            resultlines{end+1} = ['   new:  ',newlines{ii}]; %#ok<AGROW>
        end
        nshown = nshown + 1;
        if nshown >= 3
            break
        end
    end
end
r = edresult(title,description,false,resultlines);
end

function textlines = edtexlines(filename)
% The non-empty lines of a .tex file that are not comments.
textlines = strsplit(strrep(fileread(filename),sprintf('\r'),''),sprintf('\n'));
textlines = textlines(~cellfun(@isempty,strtrim(textlines)));
textlines = textlines(~strncmp(textlines,'%',1));
end

function r = edtest17(outputdirectory,filestem,~)
title = 'A change of one input parameter gives new stage files for the right stages';
description = {'1 m cube, with result recycling switched on. One reference run, and then one run per',...
    'changed parameter. For each stage file, the stored hash is compared with that of the reference',...
    'run: it should change if, and only if, the stage depends on the changed parameter.'};

stagesuffixes = {'_eddata','_Sdata','_Rdata','_paths','_tf','_ed2data','_submatrixdata','_qfinal','_tfinteq'};

% Each row: description, and which of the stage files above should get a
% new hash (1) or keep the hash (0).
changes = {
    'Receiver positions',                    [0 0 1 1 1 0 0 0 1];
    'Source position',                       [0 1 0 1 1 0 0 1 1];
    'Source amplitudes',                     [0 0 0 0 1 0 0 1 1];
    'Frequencies',                           [0 0 0 0 1 0 0 1 1];
    'Speed of sound',                        [0 0 0 0 1 0 0 1 1];
    'ngauss',                                [0 0 0 0 0 0 1 1 1];
    'difforder 10 -> 15 (GMRES)',            [0 0 0 0 0 0 0 0 0];
    'diff1method = GL',                      [0 0 0 0 1 0 0 0 0];
    'ESIEpropagationimplementation = legacy',[0 0 0 0 0 0 0 0 1];
    'savediff2result = 1',                   [0 0 0 0 0 0 0 0 1];
    };

[geoinputdata,Sinputdata,Rinputdata,controlparameters] = edsmallcube;
envdata = struct;
fh = edfilehandling(outputdirectory,filestem);
fh.suppressresultrecycling = 0;

% The reference run is computed without recycling, so that all its stage
% files are written with this filestem (with recycling, a matching
% _tfinteq file from an earlier run could be reused, and then no _qfinal
% file would be written).
fhref = fh;
fhref.suppressresultrecycling = 1;
EDmain_convex(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,fhref);
refhashes = edreadstagehashes(outputdirectory,filestem,stagesuffixes);

allok = true;
resultlines = {['Stage files: ',strjoin(stagesuffixes,' ')]};
for ii = 1:size(changes,1)
    g = geoinputdata; S = Sinputdata; R = Rinputdata; e = envdata; cp = controlparameters; f = fh;
    switch ii
        case 1,  R.coordinates(1,:) = R.coordinates(1,:) + [0.1 0 0];
        case 2,  S.coordinates = S.coordinates + [0 0.1 0];
        case 3,  S.sourceamplitudes = 2;
        case 4,  cp.frequencies = [100 250];
        case 5,  e.cair = 340;
        case 6,  cp.ngauss = 10;
        case 7,  cp.difforder = 15;
        case 8,  cp.diff1method = 'GL';
        case 9,  cp.ESIEpropagationimplementation = 'legacy';
        case 10, f.savediff2result = 1;
    end
    EDmain_convex(g,S,R,e,cp,f);
    newhashes = edreadstagehashes(outputdirectory,filestem,stagesuffixes);
    haschanged = double(~strcmp(newhashes,refhashes));
    expected = changes{ii,2};
    if isequal(haschanged,expected)
        verdict = 'OK';
    else
        verdict = ['WRONG, expected ',mat2str(expected)];
        allok = false;
    end
    resultlines{end+1} = sprintf('%-40s new hash: %s  %s',changes{ii,1},mat2str(haschanged),verdict); %#ok<AGROW>
end

r = edresult(title,description,allok,resultlines);
end

function r = edtest18(outputdirectory,filestem,~)
title = 'Reuse of the integral equation solution for new receivers';
tolreuse = 1e-12;
description = {'1 m cube, with result recycling switched on. After a first run, the receivers are moved.',...
    'The second run should reuse the _qfinal file (the integral equation is not solved again),',...
    ['and give the same result as a run without recycling, to ',num2str(tolreuse),' (relative to max abs(tftot)).']};

[geoinputdata,Sinputdata,Rinputdata,controlparameters] = edsmallcube;
fh = edfilehandling(outputdirectory,filestem);
EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,fh);   % without recycling
fh.suppressresultrecycling = 0;

Rnew = Rinputdata;
Rnew.coordinates = [-1.6 0.8 0.3; 0.3 -2.2 -0.6];
EDres2 = EDmain_convex(geoinputdata,Sinputdata,Rnew,struct,controlparameters,fh);
qfinalfile = '';
if isfield(EDres2,'inteqdiagnostics') && isfield(EDres2.inteqdiagnostics,'qfinalfile')
    qfinalfile = EDres2.inteqdiagnostics.qfinalfile;
end
reused = strncmp(qfinalfile,'reused',6);

fhnorecycling = edfilehandling(outputdirectory,[filestem,'_norecycling']);
EDres3 = EDmain_convex(geoinputdata,Sinputdata,Rnew,struct,controlparameters,fhnorecycling);
relerr = max(abs(EDres2.tftot(:)-EDres3.tftot(:)))/max(abs(EDres3.tftot(:)));

r = edresult(title,description,reused && relerr < tolreuse,...
    {['Second run, integral equation solution: ',qfinalfile],...
     ['Difference from a run without recycling: ',num2str(relerr,3)]});
end

function hashes = edreadstagehashes(outputdirectory,filestem,stagesuffixes)
% Reads the EDinputdatahash from each stage file. A missing file gives ''.
hashes = cell(1,numel(stagesuffixes));
for ii = 1:numel(stagesuffixes)
    stagefile = fullfile(outputdirectory,[filestem,stagesuffixes{ii},'.mat']);
    hashes{ii} = '';
    if exist(stagefile,'file') == 2
        filedata = load(stagefile,'EDinputdatahash');
        if isfield(filedata,'EDinputdatahash')
            hashes{ii} = filedata.EDinputdatahash;
        end
    end
end
end
