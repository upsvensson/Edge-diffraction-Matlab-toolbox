function passtest = EDverify(outputdirectory,runtest,showtext,plotdiagrams)
% EDverify runs EDtoolbox for one or several defined tests,
% and compares the results to expected results. A logfile is written.
%
% Input parameters:
%   outputdirectory (optional) The directory were the logfile will be
%                   stored. If this parameter is not given any value, no
%                   logfile will be written, but the user will get a window
%                   to specify a directory where the temporary calculation
%                   results will be stored.
%   runtest         (optional) A vector of 0 or 1, which for each position n
%                   tells if test n should be run. Default value: all 1.
%   showtext        (optional) 0 or 1; determines if the results should be printed on
%                   the screen. The results are always written to the
%                   logfile. Default value: 0.
%   plotdiagrams    (optional) 0 or 1: determines if result plots will be
%                   generated. Default value: 0.
% 
% Output parameter:
%   passtest        A vector with -1 (fail), 0 (not run), or 1 (pass) for
%   all the tests.
%
% Tests:
% 1. Response at cuboid surface, plane wave incidence, 0 Hz, no singularity
% problem. HOD test.
% 2. Field continuity across zone boundaries, single edge, 100 Hz.
% Perpendicular edge hit. Diff1 test.
% 3. Field continuity across zone boundaries, single edge, 100 Hz.
% Skewed edge hit. Diff1 test.
% 4. Field continuity across corner zone boundaries, single edge, 100 Hz.
% Diff1 test.
% 5. Field continuity for receivers close to a single edge, 100 Hz.
% Diff1 test.
% 6. Replicate a non-centered internal monopole, at 0.1 Hz.
% 7. Direct sound obscuring for a corner-on hit of an octahedron.
% 8. Direct sound obscuring for an edge-on hit of a cube.
% 
% Peter Svensson 18 April 2024 (peter.svensson@ntnu.no)
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

ntests = 8;

if nargin < 1 
    outputdirectory = [];
    savelogfile = 0;
else
    savelogfile = 1;
end
if nargin < 2
    runtest = ones(1,ntests);
end
if nargin < 3
    showtext = 0;
end
if nargin < 4
    plotdiagrams = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outputdirectory = uigetdir('', 'Select the directory for temporary calculation results');

passtest = zeros(1,ntests);

[EDversionnumber,changedate,changetime] = EDgetversion;
NN = num2str(EDversionnumber);

clockvec = clock;
datevec = date;
datevec = datevec(datevec~='-');
datetimevec = [datevec,'_',sprintf('%02d',clockvec(4)),'h', sprintf('%02d',clockvec(5)),'m',sprintf('%02d',round(clockvec(6)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a logfile and open it

if savelogfile ==1
    compstr = computer;
    compstr = lower(compstr(1:3));
    if compstr == 'mac'  
        lineending = 13;
    elseif compstr == 'sun' || compstr == 'sol'            
        lineending = 10;
    elseif compstr == 'pcw'
        lineending = [13,10];
    else
        error('ERROR: Not implemented for this computer type yet')	
    end

    mfile = mfilename('fullpath');
    [infilepath,filestem] = fileparts(mfile);
    logfilename = [infilepath,filesep,'EDverify_v',NN,'_',datetimevec,'.txt'];

    fid = fopen(logfilename,'w');
    if fid == -1
        disp('The logfile is not possible to open - check that it isn''t opened by any program!')
        return
    end
    fwrite(fid,['EDverify, EDtoolbox v. ',NN,' (last change on ',changedate,'), run on ',datetimevec,lineending],'char');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Test 1
%
% Response at surface, plane wave incidence, DC, no singularity

if runtest(1) == 1
    
    itest = 1;
    II = '1';

    if showtext > 0
        disp(' ')
        disp(['EDverify, EDtoolbox v. ',NN,' (last change on ',changedate,'), run on ',datetimevec])
        disp(' ')
        disp('*********************************************************************')
        disp(['Test ',II,': EDmain_convex, DC response at surface, plane wave incidence']);
        disp('Expected rel. errors for receiver positions 1 (in front, at surface) and 2 (behind, at surface), with')
        disp('Gauss-Legendre quadrature, and 96 points for the longest edge, and diffraction order 20, are')
        disp('1.3e-5 and 1.5e-5, respectively.')
    else
        disp(' ')
        disp(['EDverify, EDtoolbox v. ',NN,' (last change on ',changedate,'), run on ',datetimevec])
        disp(' ')
        disp('*********************************************************************')
        disp(['Test ',II,': EDmain_convex, DC response at surface, plane wave incidence']);    
    end

    mfile = mfilename('fullpath');
    [infilepath,filestem] = fileparts(mfile);

    corners = [     -0.2000   -0.4400   -0.3200
        0.2000   -0.4400   -0.3200
        0.2000    0.2000   -0.3200
       -0.2000    0.2000   -0.3200
       -0.2000   -0.4400         0
        0.2000   -0.4400         0
        0.2000    0.2000         0
       -0.2000    0.2000         0];

    planecorners = [   1     4     3     2
         5     6     7     8
         1     2     6     5
         3     4     8     7
         2     3     7     6
         1     5     8     4];

    geoinputdata = struct('corners',corners,'planecorners',planecorners);
    Sinputdata = struct('coordinates',1e8*[1 1 1]/sqrt(3));
    Rinputdata = struct('coordinates',[0 0 0.00001;0 0 -0.32001]);
    controlparameters = struct('frequencies',0);
    controlparameters.difforder = 20;
    controlparameters.docalctf = 1;
    controlparameters.docalcir = 0;
    controlparameters.docalctf_ESIEBEM = 0;
    filehandlingparameters = struct('outputdirectory',outputdirectory);
    filehandlingparameters.filestem = filestem;
    filehandlingparameters.savelogfile = 1;
    filehandlingparameters.showtext = 0;
    controlparameters.ngauss = 96;    
    controlparameters.discretizationtype = 2;
    controlparameters.Rstart = norm(Sinputdata.coordinates);
    envdata.cair = 344;

    EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);

    eval(['load ''',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tfinteq.mat'''])
    eval(['load ''',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tf.mat''',' tfdirect tfgeom tfdiff EDversionnumber'])

    tftot = tfdirect + tfgeom + tfdiff + tfinteqdiff;
    tftot = tftot*norm(Sinputdata.coordinates);

    relerr = abs(tftot-1);
    if relerr(1) < 0.14e-4 && relerr(2) < 0.16e-4
        passtest(itest) = 1;
    else
        passtest(itest) = -1;        
    end

    if showtext > 0
        disp(' ')
        disp(['Computed results are:   ',num2str(relerr(1)),' and ',num2str(relerr(2))])
        disp(' ')
        if passtest(itest) == 1
            disp(['So, verification test ',II,' was passed'])
        else
            disp(['So, verification test ',II,' was not passed. Please check the code'])
            disp('Make sure these settings were set in this script: difforder = 2, ngauss = 96, discretizationtype = 2')
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if savelogfile == 1
        fwrite(fid,[' ',lineending],'char');
        fwrite(fid,['####################################################################',lineending],'char');
        fwrite(fid,['Test ',II,': EDmain_convex, DC response at surface, plane wave incidence',lineending],'char');
        fwrite(fid,['Expected rel. errors for receiver positions 1 (in front, at surface) and 2 (behind, at surface), with',lineending],'char');
        fwrite(fid,['Gauss-Legendre quadrature, and 96 points for the longest edge, and diffraction order 20, are',lineending],'char');
        fwrite(fid,['1.3e-5 and 1.5e-5, respectively.',lineending],'char');
        fwrite(fid,[' ',lineending],'char');
        fwrite(fid,['Computed results are: ',num2str(relerr(1)),' and ',num2str(relerr(2)),lineending],'char');
        if relerr(1) < 0.14e-4 && relerr(2) < 0.16e-4
            fwrite(fid,['So, verification test ',II,' was passed',lineending],'char');
        else
            fwrite(fid,['So, verification test ',II,' was not passed. Please check the code'   ,lineending],'char');
            fwrite(fid,['Make sure these settings were set in this script: difforder = 2, ngauss = 96, discretizationtype = 2',lineending],'char');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Test 2: EDmain_convex, diff1 continuity across zone boundary, at 100 Hz
%          Perpendicular edge hit
%

if runtest(2) == 1
    itest = 2;
    II = int2str(itest);
    
    if showtext > 0
        disp(' ')
        disp('*********************************************************************')
        disp(['Test ',II,': EDmain_convex, diff1 continuity across zone boundary, at 100 Hz']);
        disp('Single edge; receivers distributed at, and very near, zone boundaries')
        disp('Perpendicular edge hit.')
        disp('Computed value at ZB should be very close to the mean value of the two')
        disp('surrounding responses (< 1e-5)')
    else
        disp(['Test ',II,': EDmain_convex, diff1 continuity across zone boundary, at 100 Hz']);        
    end
    mfile = mfilename('fullpath');
    [infilepath,filestem] = fileparts(mfile);

    corners = [     0   0   1
        0 0 -1
        -1 0 -1
        -1 0 1];

    planecorners = [4 3 2 1 ;1 2 3 4];

    geoinputdata = struct('corners',corners,'planecorners',planecorners);
    % geoinputdata.firstcornertoskip = 3;
    Sinputdata = struct('coordinates',[-1 1 0]/sqrt(2));
    phivecd = [45.05 45.01 45.001 45.0 44.999 44.99 44.95 315.05 315.01 315.001 315.0 314.999 314.99 314.95].';
    receivers = [cosd(phivecd) sind(phivecd) 0*phivecd];
    Rinputdata = struct('coordinates',receivers);
    controlparameters = struct('frequencies',100);
    controlparameters.difforder = 1;
    controlparameters.docalctf = 1;
    controlparameters.docalcir = 0;
    controlparameters.docalctf_ESIEBEM = 0;
    controlparameters.ngauss = 16;
    controlparameters.discretizationtype = 2;
    filehandlingparameters = struct('outputdirectory',outputdirectory);
    filehandlingparameters.filestem = filestem;
    filehandlingparameters.savelogfile = 1;
    filehandlingparameters.savelogfile = 1;
    filehandlingparameters.showtext = 0;
    controlparameters.Rstart = 0;

    envdata.cair = 344;

    EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);

    eval(['load ''',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tf.mat''',' tfdirect tfgeom tfdiff EDversionnumber'])

    tftot = tfdirect + tfgeom + tfdiff;

    iv1 = 8:14;
    ivmid1 = ceil(length(iv1)/2);
    tfZBdirect = tfdirect(iv1) + tfdiff(iv1) ;
    tfZBdirectmean = (   abs(tfZBdirect(ivmid1-1))   +   abs(tfZBdirect(ivmid1+1))     )/2;
    relerrZBdirect = abs( (   abs(tfZBdirectmean)-abs(tfZBdirect(ivmid1))   )/abs(tfZBdirectmean) );

    iv2 = 1:7;
    ivmid2 = ceil(length(iv2)/2);   
    tfZBspec = tfgeom(iv2) + tfdiff(iv2) ;
    tfZBspecmean = (  abs(tfZBspec(ivmid2-1))   +   abs(tfZBspec(ivmid2+1))    )/2;
    relerrZBspec = abs( ( abs(tfZBspecmean)-abs(tfZBspec(ivmid2)))/abs(tfZBspecmean) );

    if relerrZBdirect < 1e-5 & relerrZBspec < 1e-5
       passtest(itest) = 1; 
    else
        passtest(itest) = -1;        
    end

    if showtext > 0
        disp(' ')
        disp(['Computed results at the direct sound ZB differ by: ',num2str(relerrZBdirect)])
        disp(['Computed results at the specular reflection ZB differ by: ',num2str(relerrZBspec)])
        disp(' ')
        if passtest(itest) == 1
            disp(['So, verification test ',II,' was passed'])
        else
            disp(['So, verification test ',II,' was not passed. Please check the code'])   
        end
    end

    if plotdiagrams == 1
       figure
        plot(phivecd(iv1),abs(tfZBdirect),'-o',phivecd(iv1(ivmid1)),abs(tfZBdirect(ivmid1)),'*')
       xlabel('Receiver angle   [deg.]');
       ylabel('TF magnitude   [-]');
       title(['Test ',II,' Direct sound amplitude across zone boundary, perpendicular edge hit']);
       figure
       plot(phivecd(iv2),abs(tfZBspec),'-o',phivecd(iv2(ivmid2)),abs(tfZBspec(ivmid2)),'*')
       xlabel('Receiver angle   [deg.]');
       ylabel('TF magnitude   [-]');
       title(['Test ',II,' Specular reflection amplitude across zone boundary, perpendicular edge hit']);
       
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if savelogfile == 1
        fwrite(fid,[' ',lineending],'char');
        fwrite(fid,['####################################################################',lineending],'char');
        fwrite(fid,['Test ',II,': EDmain_convex, diff1 continuity across zone boundary, at 100 Hz',lineending],'char');
        fwrite(fid,['Single edge; receivers distributed at, and very near, zone boundaries',lineending],'char');
        fwrite(fid,['Perpendicular edge hit.',lineending],'char');
        fwrite(fid,['Computed value at ZB should be very close to the mean value of the two',lineending],'char');
        fwrite(fid,['surrounding responses (< 1e-5)',lineending],'char');
        fwrite(fid,[' ',lineending],'char');
        fwrite(fid,['Computed results at the direct sound ZB differ by: ',num2str(relerrZBdirect),lineending],'char');
        fwrite(fid,['Computed results at the specular reflection ZB differ by: ',num2str(relerrZBspec),lineending],'char');
        fwrite(fid,[' ',lineending],'char');
        if passtest(4) == 1
            fwrite(fid,['So, verification test ',II,' was passed',lineending],'char');
        else
            fwrite(fid,['So, verification test ',II,' was not passed. Please check the code'   ,lineending],'char');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Test 3: EDmain_convex, diff1 continuity across zone boundary, at 100 Hz
%          Skewed edge hit

if runtest(3) == 1
    itest = 3;
    II = int2str(itest);
    
    if showtext > 0
        disp(' ')
        disp('*********************************************************************')
        disp(['Test ',II,': EDmain_convex, diff1 continuity across zone boundary, at 100 Hz']);
        disp('Single edge; receivers distributed at, and very near, zone boundaries')
        disp('Skewed edge hit.')
        disp('Computed value at ZB should be very close to the mean value of the two')
        disp('surrounding responses (< 1e-5)')
    else
        disp(['Test ',II,': EDmain_convex, diff1 continuity across zone boundary, at 100 Hz']);        
    end
    mfile = mfilename('fullpath');
    [infilepath,filestem] = fileparts(mfile);

    corners = [     0   0   1
        0 0 -1
        -1 0 -1
        -1 0 1];

    planecorners = [4 3 2 1 ;1 2 3 4];

    geoinputdata = struct('corners',corners,'planecorners',planecorners);
    % geoinputdata.firstcornertoskip = 3;
    Sinputdata = struct('coordinates',[-1 1 0]/sqrt(2));
    phivecd = [45.05 45.01 45.001 45.0 44.999 44.99 44.95 315.05 315.01 315.001 315.0 314.999 314.99 314.95].';
    receivers = [cosd(phivecd) sind(phivecd) 1*ones(size(phivecd))];
    Rinputdata = struct('coordinates',receivers);
    controlparameters = struct('frequencies',100);
    controlparameters.difforder = 1;
    controlparameters.docalctf = 1;
    controlparameters.docalcir = 0;
    controlparameters.docalctf_ESIEBEM = 0;
    controlparameters.ngauss = 16;
    controlparameters.discretizationtype = 2;
    filehandlingparameters = struct('outputdirectory',outputdirectory);
    filehandlingparameters.filestem = filestem;
    filehandlingparameters.savelogfile = 1;
    filehandlingparameters.showtext = 0;
    controlparameters.Rstart = 0;

    envdata.cair = 344;

    EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);

    eval(['load ''',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tf.mat''',' tfdirect tfgeom tfdiff EDversionnumber'])

    tftot = tfdirect + tfgeom + tfdiff;

    iv1 = 8:14;
    ivmid1 = ceil(length(iv1)/2);
    tfZBdirect = tfdirect(iv1) + tfdiff(iv1) ;
    tfZBdirectmean = (   abs(tfZBdirect(ivmid1-1))   +   abs(tfZBdirect(ivmid1+1))     )/2;
    relerrZBdirect = abs( (   abs(tfZBdirectmean)-abs(tfZBdirect(ivmid1))   )/abs(tfZBdirectmean) );

    iv2 = 1:7;
    ivmid2 = ceil(length(iv2)/2);   
    tfZBspec = tfgeom(iv2) + tfdiff(iv2) ;
    tfZBspecmean = (  abs(tfZBspec(ivmid2-1))   +   abs(tfZBspec(ivmid2+1))    )/2;
    relerrZBspec = abs( ( abs(tfZBspecmean)-abs(tfZBspec(ivmid2)))/abs(tfZBspecmean) );

    if relerrZBdirect < 1e-5 & relerrZBspec < 1e-5
       passtest(itest) = 1; 
    else
        passtest(itest) = -1;        
    end

    if showtext > 0
        disp(' ')
        disp(['Computed results at the direct sound ZB differ by: ',num2str(relerrZBdirect)])
        disp(['Computed results at the specular reflection ZB differ by: ',num2str(relerrZBspec)])
        disp(' ')
        if passtest(itest) == 1
            disp(['So, verification test ',II,' was passed'])
        else
            disp(['So, verification test ',II,' was not passed. Please check the code'])   
        end
    end

    if plotdiagrams == 1
       figure
        plot(phivecd(iv1),abs(tfZBdirect),'-o',phivecd(iv1(ivmid1)),abs(tfZBdirect(ivmid1)),'*')
       xlabel('Receiver angle   [deg.]');
       ylabel('TF magnitude   [-]');
       title(['Test ',II,' Direct sound amplitude across zone boundary, skewed edge hit']);
       figure
       plot(phivecd(iv2),abs(tfZBspec),'-o',phivecd(iv2(ivmid2)),abs(tfZBspec(ivmid2)),'*')
       xlabel('Receiver angle   [deg.]');
       ylabel('TF magnitude   [-]');
       title(['Test ',II,' Specular reflection amplitude across zone boundary, skewed edge hit']);
       
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if savelogfile == 1
        fwrite(fid,[' ',lineending],'char');
        fwrite(fid,['####################################################################',lineending],'char');
        fwrite(fid,['Test ',II,': EDmain_convex, diff1 continuity across zone boundary, at 100 Hz',lineending],'char');
        fwrite(fid,['Single edge; receivers distributed at, and very near, zone boundaries',lineending],'char');
        fwrite(fid,['Skewed edge hit.',lineending],'char');
        fwrite(fid,['Computed value at ZB should be very close to the mean value of the two',lineending],'char');
        fwrite(fid,['surrounding responses (< 1e-5)',lineending],'char');
        fwrite(fid,[' ',lineending],'char');
        fwrite(fid,['Computed results at the direct sound ZB differ by: ',num2str(relerrZBdirect),lineending],'char');
        fwrite(fid,['Computed results at the specular reflection ZB differ by: ',num2str(relerrZBspec),lineending],'char');
        fwrite(fid,[' ',lineending],'char');
        if passtest(4) == 1
            fwrite(fid,['So, verification test ',II,' was passed',lineending],'char');
        else
            fwrite(fid,['So, verification test ',II,' was not passed. Please check the code'   ,lineending],'char');
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Test 4: EDmain_convex, diff1 continuity across corner zone boundary, at 100 Hz
%

if runtest(4) == 1
    itest = 4;
    II = int2str(itest);
    
    if showtext > 0
        disp(' ')
        disp('*********************************************************************')
        disp(['Test ',II,': EDmain_convex, diff1 continuity across corner zone boundary, at 100 Hz']);
        disp('Two edges meat at 90 deg.; several receivers exactly at, and very near zone boundaries')
        disp('Computed value at ZB should be very close to the mean value of the two')
        disp('surrounding responses (< 1e-3)')
    else
        disp(['Test ',II,': EDmain_convex, diff1 continuity across corner zone boundary, at 100 Hz']);    
    end

    mfile = mfilename('fullpath');
    [infilepath,filestem] = fileparts(mfile);

    corners = [     0   0   -10
        0 0 1
        -10 0 1
        -10 0 -10];

    planecorners = [4 3 2 1 ;1 2 3 4];

    geoinputdata = struct('corners',corners,'planecorners',planecorners);
    geoinputdata.firstcornertoskip = 4;
    Sinputdata = struct('coordinates',[-1 -1 0]);
    receivers =[0.99 -1 1.98;
                0.995 -1 1.99; 
                1 -1 2; 
                1.005 -1 2.01;
                1.01 -1 2.02;
                0.99 1 1.98
                0.995 1 1.99;
                1 1 2; 
                1.005 1 2.01
                1.01 1 2.02];
    Rinputdata = struct('coordinates',receivers);
    controlparameters = struct('frequencies',100);
    controlparameters.difforder = 4;
    controlparameters.docalctf = 1;
    controlparameters.docalcir = 0;
    controlparameters.docalctf_ESIEBEM = 0;
    controlparameters.discretizationtype = 2;
    controlparameters.ngauss = 32;
    filehandlingparameters = struct('outputdirectory',outputdirectory);
    filehandlingparameters.filestem = filestem;
    filehandlingparameters.savelogfile = 1;
    filehandlingparameters.showtext = 0;
    controlparameters.Rstart = 0;

    envdata.cair = 344;

    EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);

    eval(['load ''',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tf.mat''',' tfdirect tfgeom tfdiff EDversionnumber'])
    eval(['load ''',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tfinteq.mat'''])

    tftot = tfdirect + tfgeom + tfdiff;

    iv = 6:10;
    ivmid = ceil(length(iv)/2);
    tfZBdirect = tfdirect(iv) + tfdiff(iv) + tfinteqdiff(iv);
    tfZBdirectmean = (   abs(tfZBdirect(ivmid-1))   +   abs(tfZBdirect(ivmid+1))     )/2;
    relerrZBdirect = abs( (   abs(tfZBdirectmean)-abs(tfZBdirect(ivmid))   )/abs(tfZBdirectmean) );

    iv = 1:5;
    ivmid = ceil(length(iv)/2);   
    tfZBspec = tfgeom(iv) + tfdiff(iv) + tfinteqdiff(iv);
    tfZBspecmean = (  abs(tfZBspec(ivmid-1))   +   abs(tfZBspec(ivmid+1))    )/2;
    relerrZBspec = abs( ( abs(tfZBspecmean)-abs(tfZBspec(ivmid)))/abs(tfZBspecmean) );

    if relerrZBdirect < 1e-3 & relerrZBspec < 1e-3
       passtest(itest) = 1; 
    else
        passtest(itest) = -1;        
    end

    if showtext > 0
         disp(' ')
        disp(['Computed results at the direct sound corner ZB differ by: ',num2str(relerrZBdirect)])
        disp(['Computed results at the specular reflection corner ZB differ by: ',num2str(relerrZBspec)])
        disp(['   '])
        if passtest(itest) == 1
            disp(['So, verification test ',II,' was passed'])
        else
            disp(['So, verification test ',II,' was not passed. Please check the code'])   
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if savelogfile == 1
        fwrite(fid,[' ',lineending],'char');
        fwrite(fid,['####################################################################',lineending],'char');
        fwrite(fid,['Test ',II,': EDmain_convex, diff1 continuity across corner zone boundary, at 100 Hz',lineending],'char');
        fwrite(fid,['Two edges meat at 90 deg.; several receivers exactly at, and very near zone boundaries',lineending],'char');
        fwrite(fid,['Computed value at ZB should be very close to the mean value of the two',lineending],'char');
        fwrite(fid,['surrounding responses (< 1e-3)',lineending],'char');
        fwrite(fid,[' ',lineending],'char');
        fwrite(fid,['Computed results at the direct sound corner ZB differ by: ',num2str(relerrZBdirect),lineending],'char');
        fwrite(fid,['Computed results at the specular reflection corner ZB differ by: ',num2str(relerrZBspec),lineending],'char');
        fwrite(fid,[' ',lineending],'char');
        if relerrZBdirect < 1e-3 & relerrZBspec < 1e-3
            fwrite(fid,['So, verification test ',II,' was passed',lineending],'char');
        else
            fwrite(fid,['So, verification test ',II,' was not passed. Please check the code'   ,lineending],'char');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Test 5: EDmain_convex, diff1 continuity across edge, at 100 Hz
%

if runtest(5) == 1
    itest = 5;
    II = int2str(itest);
    
    if showtext > 0
        disp(' ')
        disp('*********************************************************************')
        disp(['Test ',II,': EDmain_convex, diff1 continuity across edge, for very low frequencies']);
        disp('Single edge; receivers very near the edge, in front and in back')
        disp('Far source, perp. incidence and skewed incidence.')
        disp('Computed values in front and in back should be within 1.5e-3 for the three frequencies')
    else
        disp(['Test ',II,': EDmain_convex, diff1 continuity across edge, for very low frequencies']);    
    end

    mfile = mfilename('fullpath');
    [infilepath,filestem] = fileparts(mfile);

    corners = [     0   0   -100
        0 0 100
        -10 0 100
        -10 0 -100];

    planecorners = [4 3 2 1 ;1 2 3 4];

    geoinputdata = struct('corners',corners,'planecorners',planecorners);
    geoinputdata.firstcornertoskip = 3;
    soudist = 1e3;
    sources = [[-1 -1 0]*soudist/sqrt(2);[-1 -1 1]*soudist/sqrt(3)];
    Sinputdata = struct('coordinates',sources);
    receivers = [-0.000001 -0.00000001 0;-0.000001 0.00000001 0];
    Rinputdata = struct('coordinates',receivers);
    controlparameters = struct('frequencies',[0.1 ]);
    controlparameters.difforder = 1;
    controlparameters.docalctf = 1;
    controlparameters.docalcir = 0;
    controlparameters.docalctf_ESIEBEM = 0;
    controlparameters.discretizationtype = 2;
    controlparameters.ngauss = 26;
    filehandlingparameters = struct('outputdirectory',outputdirectory);
    filehandlingparameters.filestem = filestem;
    filehandlingparameters.savelogfile = 1;
    filehandlingparameters.showtext = 0;
    controlparameters.Rstart = soudist;

    envdata.cair = 344;

    EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);

    eval(['load ''',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tf.mat''',' tfdirect tfgeom tfdiff EDversionnumber'])

    tftot = soudist*(tfdirect + tfgeom + tfdiff);
    
    tftot = [squeeze(tftot(:,:,1)) squeeze(tftot(:,:,2))];
    
    pressurediff = [abs(tftot(:,1) - tftot(:,2));abs(tftot(:,3) - tftot(:,4))];
    
    if max(pressurediff) < 2.3e-3 
       passtest(itest) = 1; 
    else
        passtest(itest) = -1;        
    end

    if showtext > 0
         disp(' ')
        disp(['Computed results differ by max: ',num2str(max(pressurediff))])
        disp(['   '])
        if passtest(itest) == 1
            disp(['So, verification test ',II,' was passed'])
        else
            disp(['So, verification test ',II,' was not passed. Please check the code'])   
        end
    end

    if plotdiagrams == 1
       figure
       semilogx(controlparameters.frequencies,abs(tftot),'-o')       
        grid
       xlabel('Frequency   [Hz]');
       ylabel('TF magnitude   [-]');
       title(['Test ',II,' Sound pressure amplitude very near edge (',num2str(norm(receivers(1,:))),' m from the edge)']);  
       legend('Frontal side, perp.incidence','Shadow side, perp.incidence','Frontal side, skewed incidence','Shadow side, skewed incidence')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if savelogfile == 1
        fwrite(fid,[' ',lineending],'char');
        fwrite(fid,['####################################################################',lineending],'char');
        fwrite(fid,['Test ',II,': EDmain_convex, diff1 continuity across edge, for very low frequencies',lineending],'char');
        fwrite(fid,['Single edge; receivers very near the edge, in front and in back',lineending],'char');
        fwrite(fid,['Far source, perp. incidence and skewed incidence.',lineending],'char');
        fwrite(fid,['Computed values in front and in back should be within 1.5e-3 for the three frequencies',lineending],'char');
        fwrite(fid,['Computed results differ by max: ',num2str(max(pressurediff)),lineending],'char');
        fwrite(fid,[' ',lineending],'char');
        if passtest(itest) == 1
            fwrite(fid,['So, verification test ',II,' was passed',lineending],'char');
        else
            fwrite(fid,['So, verification test ',II,' was not passed. Please check the code'   ,lineending],'char');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Test 6: EDmain_convex, replicate internal monopole, at 0.1 Hz

if runtest(6) == 1
    itest = 6;
    II = int2str(itest);
    
    if showtext > 0
        disp(' ')
        disp('*********************************************************************')
        disp(['Test ',II,': EDmain_convex, replicate internal monopole, at 0.1 Hz']);
        disp('Cube w internal, non-centered monopole')
        disp('Entire cube surface gets sources')
        disp('Radiated field should be within [0.996,1.008] around a circle of receivers')
    else
        disp(['Test ',II,': EDmain_convex, replicate internal monopole, at 0.1 Hz']);    
    end

    mfile = mfilename('fullpath');
    [infilepath,filestem] = fileparts(mfile);

    corners = [     -0.5000   -0.500   -0.500
        0.5000   -0.500   -0.500
        0.5000    0.5000   -0.500
       -0.5000    0.5000   -0.500
       -0.5000   -0.500         0.5
        0.5000   -0.500         0.5
        0.5000    0.5000         0.5
       -0.5000    0.5000         0.5];
    planecorners = [   1     4     3     2
         5     6     7     8
         1     2     6     5
         3     4     8     7
         2     3     7     6
         1     5     8     4];
    geoinputdata = struct('corners',corners);
    geoinputdata.planecorners = planecorners;
    planedata = EDreadgeomatrices(geoinputdata.corners,geoinputdata.planecorners);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Controlparameters

    controlparameters = struct('frequencies',0.1);
    nfrequencies = length(controlparameters.frequencies);
    controlparameters.ngauss = 32;
    controlparameters.difforder = 15;
    controlparameters.docalctf = 1;
    controlparameters.docalcir = 0;
    controlparameters.docalctf_ESIEBEM = 0;
    controlparameters.discretizationtype = 2;

    envdata = struct('cair',344);
    envdata.rhoair = 1.21;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sources

    patchpoints = [];
    patchcorners = [];
    patchareas = [];
    patchnvecs = [];
    n1 = 6;
    n2 = 6;
    planenvecs = planedata.planeeqs(:,1:3);
    for ii = 1:6
       [patchcornercoords,patchmidcoords,areas] = EDdivrect(corners(planecorners(ii,1),:),corners(planecorners(ii,2),:),corners(planecorners(ii,3),:),corners(planecorners(ii,4),:),n1,n2);     
       nvec = planenvecs(ii,:);
       patchmidcoords = patchmidcoords + 0.0001*nvec(ones(n1*n2,1),:);
       patchpoints = [patchpoints;patchmidcoords];
       patchcorners = [patchcorners;patchcornercoords];
       patchareas = [patchareas;areas];
       patchnvecs = [patchnvecs;nvec(ones(n1*n2,1),:)];
    end

    k = 2*pi*controlparameters.frequencies/envdata.cair;

    internalmonopole = [0.3 0.3 0.3];
    quadweights = [64 25 25 25 25 40 40 40 40]/81/4;
    quadshift = sqrt(3)/sqrt(5);
    umono_normal_patch = 0;

    % Gauss point 1 = center point
    patchquadraturepoints = patchpoints;
    [distances,cosfi] = EDcalccosfi(internalmonopole,patchquadraturepoints,patchnvecs);
    pmono_patch = exp(-1i*k*distances)./distances;
    umono_radial_patch = pmono_patch/envdata.cair/envdata.rhoair.*(1 + 1./(1i*k*distances));
    umono_normal_patch =  quadweights(1)*umono_radial_patch.*cosfi;

    % Gauss points 2-5 = from center point, towards corners
    for ii = 1:4
        patchquadraturepoints = patchpoints + quadshift*(patchcorners(:,[1:3] + 3*(ii-1)) - patchpoints);
        [distances,cosfi] = EDcalccosfi(internalmonopole,patchquadraturepoints,patchnvecs);
        pmono_patch = exp(-1i*k*distances)./distances;
        umono_radial_patch = pmono_patch/envdata.cair/envdata.rhoair.*(1 + 1./(1i*k*distances));
        umono_normal_patch = umono_normal_patch + quadweights(ii+1)*umono_radial_patch.*cosfi;    
    end

    % Gauss points 6-9 = from center point, towards midpoint
    for ii = 1:4
        if ii < 4
            patchedgemidpoint = 0.5*( patchcorners(:,[1:3] + 3*(ii-1)) + patchcorners(:,[1:3] + 3*(ii)) );
        else
            patchedgemidpoint = 0.5*( patchcorners(:,10:12) + patchcorners(:,1:3) );        
        end
        patchquadraturepoints = patchpoints + quadshift*(patchedgemidpoint - patchpoints);
        [distances,cosfi] = EDcalccosfi(internalmonopole,patchquadraturepoints,patchnvecs);
        pmono_patch = exp(-1i*k*distances)./distances;
        umono_radial_patch = pmono_patch/envdata.cair/envdata.rhoair.*(1 + 1./(1i*k*distances));
        umono_normal_patch = umono_normal_patch + quadweights(ii+5)*umono_radial_patch.*cosfi;    
    end

    Sinputdata = struct('coordinates',patchpoints);
    nsources = size(Sinputdata.coordinates,1);
    Sinputdata.sourceamplitudes = patchareas.*umono_normal_patch...
         .*1i*2*pi*controlparameters.frequencies*envdata.rhoair/4/pi;
    Sinputdata.doaddsources = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Receivers

    rdist = 3;
    phivec = [0:0.1:2*pi].';
    receivers = rdist*[cos(phivec) sin(phivec) 0*phivec];
    nreceivers = size(receivers,1);
    receivers = receivers + internalmonopole(ones(nreceivers,1),:);

    Rinputdata = struct('coordinates',receivers);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Filehandling

    filehandlingparameters = struct('outputdirectory',outputdirectory);
    filehandlingparameters.filestem = filestem;
    filehandlingparameters.showtext = 0;
    filehandlingparameters.savelogfile = 1;

    EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load and present the results

    eval(['load ''',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tfinteq.mat'''])
    eval(['load ''',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tf.mat'''])

    tftot = tfinteqdiff + tfdiff + tfdirect + tfgeom;

    maxval = rdist*max(abs(tftot));
    minval = rdist*min(abs(tftot));

    if maxval < 1.008 && minval > 0.996
       passtest(itest) = 1; 
    else
       passtest(itest) = -1;         
    end
    
    if plotdiagrams == 1
       figure

        plot(phivec*180/pi,rdist*abs([tftot.' ]),'-o')
        xlabel(['Receiver angle   [deg.]'])
        ylabel(['Sound pressure amp., re. free-field   [-]'])       
        grid
    end
    
    if showtext > 0
         disp(' ')
        disp(['Computed results are within [',num2str(minval),',',num2str(maxval)])
        disp(['   '])
        if passtest(itest) == 1
            disp(['So, verification test ',II,' was passed'])
        else
            disp(['So, verification test ',II,' was not passed. Please check the code'])   
        end
    end
    
    if savelogfile == 1
        fwrite(fid,[' ',lineending],'char');
        fwrite(fid,['####################################################################',lineending],'char');
        fwrite(fid,['Test ',II,': EDmain_convex, replicate internal monopole, at 0.1 Hz.',lineending],'char');
        fwrite(fid,['Cube w internal, non-centered monopole.',lineending],'char');
        fwrite(fid,['Entire cube surface gets sources.',lineending],'char');
        fwrite(fid,['Radiated field should be within [0.996,1.008] around a circle of receivers',lineending],'char');
        fwrite(fid,['Computed results are within [',num2str(minval),',',num2str(maxval),']',lineending],'char');
        fwrite(fid,[' ',lineending],'char');
        if passtest(itest) == 1
            fwrite(fid,['So, verification test ',II,' was passed',lineending],'char');
        else
            fwrite(fid,['So, verification test ',II,' was not passed. Please check the code'   ,lineending],'char');
        end
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Test 7 EDmain_convex, make sure the direct sound is obscured for a
%  corner-on hit, for an octahedron
%

if runtest(7) == 1
    itest = 7;
    II = int2str(itest);
    
    if showtext > 0
        disp(' ')
        disp('*********************************************************************')
        disp(['Test ',II,': EDmain_convex, make sure the direct sound is obscured for a']);
        disp('corner-on hit, for an octahedron')
    else
        disp(['Test ',II,': EDmain_convex, direct sound obscuring for corner-on hit, octahedron']);    
    end

    mfile = mfilename('fullpath');
    [infilepath,filestem] = fileparts(mfile);

    corners = [
       0.707106781 0 0
        0	0.707106781	0
        -0.707106781	0	0
        0	-0.707106781	0
        0	0	0.707106781
        0	0	-0.707106781];

    planecorners = [
    1 2 5
    2 3 5
    3 4 5
    4 1 5
    2 1 6
    3 2 6
    4 3 6
    1 4 6 ];

    geoinputdata = struct('corners',corners,'planecorners',planecorners);
    geoinputdata.firstcornertoskip = 1e6;
    Sinputdata = struct('coordinates',[2 0 0]);
    receivers = [-2 0 0];
    Rinputdata = struct('coordinates',receivers);
    controlparameters = struct('frequencies',100);
    controlparameters.difforder = 0;
    controlparameters.docalctf = 1;
    controlparameters.docalcir = 0;
    controlparameters.docalctf_ESIEBEM = 0;
    controlparameters.discretizationtype = 2;
    controlparameters.ngauss = 16;
    
    filehandlingparameters = struct('outputdirectory',outputdirectory);
    filehandlingparameters.filestem = filestem;
    filehandlingparameters.savelogfile = 1;
    filehandlingparameters.showtext = 0;
    controlparameters.Rstart = 0;

    envdata.cair = 344;

    EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);

    eval(['load ''',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tf.mat''',' tfdirect'])
%     eval(['load ',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tfinteq.mat'])

    if abs(tfdirect) == 0
       passtest(itest) = 1; 
    else
        passtest(itest) = -1;        
    end

    if showtext > 0
         disp(' ')
        disp(['The direct sound should be zero. The computed direct sound had the value ',num2str(abs(tfdirect))])
        disp(['   '])
        if passtest(itest) == 1
            disp(['So, verification test ',II,' was passed'])
        else
            disp(['So, verification test ',II,' was not passed. Please check the code'])   
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if savelogfile == 1
        fwrite(fid,[' ',lineending],'char');
        fwrite(fid,['####################################################################',lineending],'char');
        fwrite(fid,['Test ',II,': EDmain_convex, make sure the direct sound is obscured for a',lineending],'char');
        fwrite(fid,['corner-on hit, for an octahedron',lineending],'char');
        fwrite(fid,[' ',lineending],'char');
        fwrite(fid,['The direct sound should be zero. The computed direct sound had the value ',num2str(abs(tfdirect)),lineending],'char');
        fwrite(fid,[' ',lineending],'char');
        if passtest(itest) == 1
            fwrite(fid,['So, verification test ',II,' was passed',lineending],'char');
        else
            fwrite(fid,['So, verification test ',II,' was not passed. Please check the code'   ,lineending],'char');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Test 8 EDmain_convex, make sure the direct sound is obscured for an
%  edge-on hit, for a cube
%

if runtest(8) == 1
    itest = 8;
    II = int2str(itest);
    
    if showtext > 0
        disp(' ')
        disp('*********************************************************************')
        disp(['Test ',II,': EDmain_convex, make sure the direct sound is obscured for an']);
        disp('edge-on hit, for a cube')
    else
        disp(['Test ',II,': EDmain_convex, direct sound obscuring for edge-on hit, cube']);    
    end


    mfile = mfilename('fullpath');
    [infilepath,filestem] = fileparts(mfile);

    corners = [     -0.5000   -0.500   -0.500
        0.5000   -0.50   -0.500
        0.5000    0.5000   -0.500
       -0.5000    0.5000   -0.500
       -0.5000   -0.500         0.5
        0.5000   -0.500         0.5
        0.5000    0.5000         0.5
       -0.5000    0.5000         0.5];

    planecorners = [   1     4     3     2
         5     6     7     8
         1     2     6     5
         3     4     8     7
         2     3     7     6
         1     5     8     4];

    geoinputdata = struct('corners',corners,'planecorners',planecorners);
    sources = [1 1 0];
    Sinputdata = struct('coordinates',sources);
    receivers = [-1 -1 0];
    Rinputdata = struct('coordinates',receivers);
    controlparameters = struct('frequencies',100);
    controlparameters.difforder = 0;
    controlparameters.docalctf = 1;
    controlparameters.docalcir = 0;
    controlparameters.docalctf_ESIEBEM = 0;
    controlparameters.discretizationtype = 2;
    controlparameters.ngauss = 16;
    
    filehandlingparameters = struct('outputdirectory',outputdirectory);
    filehandlingparameters.filestem = filestem;
    filehandlingparameters.savelogfile = 1;
    filehandlingparameters.showtext = 0;
    controlparameters.Rstart = 0;

    EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load and present the results

    eval(['load ''',outputdirectory,filesep,filehandlingparameters.filestem,'_tf.mat''',' tfdirect'])

    tfdirect
    
   if abs(tfdirect) == 0
       passtest(itest) = 1; 
    else
        passtest(itest) = -1;        
    end

    if showtext > 0
         disp(' ')
        disp(['The direct sound should be zero. The computed direct sound had the value ',num2str(abs(tfdirect))])
        disp(['   '])
        if passtest(itest) == 1
            disp(['So, verification test ',II,' was passed'])
        else
            disp(['So, verification test ',II,' was not passed. Please check the code'])   
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if savelogfile == 1
        fwrite(fid,[' ',lineending],'char');
        fwrite(fid,['####################################################################',lineending],'char');
        fwrite(fid,['Test ',II,': EDmain_convex, make sure the direct sound is obscured for an',lineending],'char');
        fwrite(fid,['edge-on hit, for a cube',lineending],'char');
        fwrite(fid,[' ',lineending],'char');
        fwrite(fid,['The direct sound should be zero. The computed direct sound had the value ',num2str(abs(tfdirect)),lineending],'char');
        fwrite(fid,[' ',lineending],'char');
        if passtest(itest) == 1
            fwrite(fid,['So, verification test ',II,' was passed',lineending],'char');
        else
            fwrite(fid,['So, verification test ',II,' was not passed. Please check the code'   ,lineending],'char');
        end
    end

end

if savelogfile == 1
    fclose(fid);
end

