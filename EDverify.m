function passtest = EDverify(runtest,showtext,plotdiagrams)
% EDverify runs EDtoolbox for one or several defined tests,
% and compares the results to expected results. A logfile is written.
%
% Input parameters:
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
% 1. 
% 
% Peter Svensson 16 Jan. 2018 (peter.svensson@ntnu.no)
% 
% passtest = EDverify(runtest,showtext,plotdiagrams);

% 15 Jan. 2018 First version
% 15 Jan. 2018 Cleaned up the EDversion/EDversionnumber confusion. Added a
% test 4, with fine-grain ZB check.
% 16 Jan. 2018 Changed the function call to EDgetversion instead of
% EDversion.
% 16 Jan. 2018 Changed EDverify into a function rather than a script.
%              Added the plotdiagrams parameter.

ntests = 4;

if nargin < 1
    runtest = ones(1,ntests);
end
if nargin < 2
    showtext = 0;
end
if nargin < 3
    plotdiagrams = 0;
end

passtest = zeros(1,ntests);

[EDversionnumber,changedate,changetime] = EDgetversion;
NN = num2str(EDversionnumber);

clockvec = clock;
datevec = date;
datevec = datevec(datevec~='-');
datetimevec = [datevec,'_',sprintf('%02d',clockvec(4)),'h', sprintf('%02d',clockvec(5)),'m',sprintf('%02d',round(clockvec(6)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a logfile and open it

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
logfilename = [infilepath,filesep,'results',filesep,'EDverify_v',NN,'_',datetimevec,'.txt'];

fid = fopen(logfilename,'w');
if fid == -1
    disp('The logfile is not possible to open - check that it isn''t opened by any program!')
    return
end
fwrite(fid,['EDverify, EDtoolbox v. ',NN,' (last change on ',changedate,'), run on ',datetimevec,lineending],'char');

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
        disp(['Test ',II,': EDmain_convexESIE, DC response at surface, plane wave incidence']);
        disp('Expected rel. errors for receiver positions 1 (in front, at surface) and 2 (behind, at surface), with')
        disp('Gauss-Legendre quadrature, and 96 points for the longest edge, and diffraction order 20, are')
        disp('1.3e-5 and 1.5e-5, respectively.')
    else
        disp(' ')
        disp(['EDverify, EDtoolbox v. ',NN,' (last change on ',changedate,'), run on ',datetimevec])
        disp(' ')
        disp('*********************************************************************')
        disp(['Test ',II,': EDmain_convexESIE, DC response at surface, plane wave incidence']);    
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

    geofiledata = struct('corners',corners,'planecorners',planecorners);
    Sindata = struct('coordinates',1e8*[1 1 1]/sqrt(3));
    Rindata = struct('coordinates',[0 0 0.00001;0 0 -0.32001]);
    controlparameters = struct('frequencies',0);
    controlparameters.difforder = 20;
    filehandlingparameters = struct('outputdirectory',infilepath);
    filehandlingparameters.filestem = filestem;
    filehandlingparameters.savelogfile = 1;
    filehandlingparameters.showtext = 0;
    controlparameters.ngauss = 96;    
    controlparameters.discretizationtype = 2;
    controlparameters.Rstart = norm(Sindata.coordinates);
    envdata.cair = 344;

    EDmain_convexESIE(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);

    eval(['load ',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tfinteq.mat'])
    eval(['load ',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tf.mat tfdirect tfgeom tfdiff EDversionnumber'])

    tftot = tfdirect + tfgeom + tfdiff + tfinteqdiff;
    tftot = tftot*norm(Sindata.coordinates);

    relerr = abs(tftot-1);
    if relerr(1) < 0.14e-4 && relerr(2) < 0.16e-4
        passtest(itest) = 1;
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

    fwrite(fid,[' ',lineending],'char');
    fwrite(fid,['####################################################################',lineending],'char');
    fwrite(fid,['Test ',II,': EDmain_convexESIE, DC response at surface, plane wave incidence',lineending],'char');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Test 2: EDmain_convexESIE, diff1 continuity across zone boundary, at 100 Hz
%

if runtest(2) == 1
    itest = 2;
    II = int2str(itest);
    
    if showtext > 0
        disp(' ')
        disp('*********************************************************************')
        disp(['Test ',II,': EDmain_convexESIE, diff1 continuity across zone boundary, at 100 Hz']);
        disp('Single edge; three receivers exactly at, and very near zone boundaries')
        disp('Computed value at ZB should be very close to the mean value of the two')
        disp('surrounding responses (< 1e-5)')
    else
        disp(['Test ',II,': EDmain_convexESIE, diff1 continuity across zone boundary, at 100 Hz']);        
    end

    mfile = mfilename('fullpath');
    [infilepath,filestem] = fileparts(mfile);

    corners = [     0   0   1
        0 0 -1
        -1 0 -1
        -1 0 1];

    planecorners = [4 3 2 1 ;1 2 3 4];

    geofiledata = struct('corners',corners,'planecorners',planecorners);
%     geofiledata.firstcornertoskip = 3;
    Sindata = struct('coordinates',[-1 1 0]/sqrt(2));
    phivecd = [45.1 45.0 44.9 315.1 315.0 314.9].';
    receivers = [cosd(phivecd) sind(phivecd) 0*phivecd];
    Rindata = struct('coordinates',receivers);
    controlparameters = struct('frequencies',100);
    controlparameters.difforder = 1;
    filehandlingparameters = struct('outputdirectory',infilepath);
    filehandlingparameters.filestem = filestem;
    filehandlingparameters.savelogfile = 1;
    filehandlingparameters.savelogfile = 1;
    filehandlingparameters.showtext = 0;
    controlparameters.Rstart = 0;

    envdata.cair = 344;

    EDmain_convexESIE(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);

    eval(['load ',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tf.mat tfdirect tfgeom tfdiff EDversionnumber'])

    tftot = tfdirect + tfgeom + tfdiff;

    iv1 = 4:6;
    ivmid1 = ceil(length(iv1)/2);    
    tfZBdirect = tfdirect(iv1) + tfdiff(iv1);
    tfZBdirectmean = (   abs(tfZBdirect(ivmid1-1))   +   abs(tfZBdirect(ivmid1+1))     )/2;
    relerrZBdirect = abs( (abs(tfZBdirectmean)-abs(tfZBdirect(ivmid1)))/abs(tfZBdirectmean) );

    
    iv2 = 1:3;
    ivmid2 = ceil(length(iv2)/2);    
    tfZBspec = tfgeom(iv2) + tfdiff(iv2);
    tfZBspecmean = (  abs(tfZBspec(ivmid2-1))   +   abs(tfZBspec(ivmid2+1))    )/2;
    relerrZBspec = abs( ( abs(tfZBspecmean)-abs(tfZBspec(ivmid2)))/abs(tfZBspecmean) );

    if relerrZBdirect < 1e-5 & relerrZBspec < 1e-5
       passtest(itest) = 1; 
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
       title(['Test ',II,' Direct sound amplitude across zone boundary']);
       figure
       plot(phivecd(iv2),abs(tfZBspec),'-o',phivecd(iv2(ivmid2)),abs(tfZBspec(ivmid2)),'*')
       xlabel('Receiver angle   [deg.]');
       ylabel('TF magnitude   [-]');
       title(['Test ',II,' Specular reflection amplitude across zone boundary']);
       
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fwrite(fid,[' ',lineending],'char');
    fwrite(fid,['####################################################################',lineending],'char');
    fwrite(fid,['Test ',II,': EDmain_convexESIE, diff1 continuity across zone boundary, at 100 Hz',lineending],'char');
    fwrite(fid,['Single edge; three receivers exactly at, and very near zone boundaries',lineending],'char');
    fwrite(fid,['Computed value at ZB should be very close to the mean value of the two',lineending],'char');
    fwrite(fid,['surrounding responses (< 1e-5)',lineending],'char');
    fwrite(fid,[' ',lineending],'char');
    fwrite(fid,['Computed results at the direct sound ZB differ by: ',num2str(relerrZBdirect),lineending],'char');
    fwrite(fid,['Computed results at the specular reflection ZB differ by: ',num2str(relerrZBspec),lineending],'char');
    fwrite(fid,[' ',lineending],'char');
    if relerrZBdirect < 1e-5 & relerrZBspec < 1e-5
        fwrite(fid,['So, verification test ',II,' was passed',lineending],'char');
    else
        fwrite(fid,['So, verification test ',II,' was not passed. Please check the code'   ,lineending],'char');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Test 3: EDmain_convexESIE, diff1 continuity across corner zone boundary, at 100 Hz
%

if runtest(3) == 1
    itest = 3;
    II = int2str(itest);
    
    if showtext > 0
        disp(' ')
        disp('*********************************************************************')
        disp(['Test ',II,': EDmain_convexESIE, diff1 continuity across corner zone boundary, at 100 Hz']);
        disp('Two edges meat at 90 deg.; three receivers exactly at, and very near zone boundaries')
        disp('Computed value at ZB should be very close to the mean value of the two')
        disp('surrounding responses (< 1e-3)')
    else
        disp(['Test ',II,': EDmain_convexESIE, diff1 continuity across corner zone boundary, at 100 Hz']);    
    end

    mfile = mfilename('fullpath');
    [infilepath,filestem] = fileparts(mfile);

    corners = [     0   0   -10
        0 0 1
        -10 0 1
        -10 0 -10];

    planecorners = [4 3 2 1 ;1 2 3 4];

    geofiledata = struct('corners',corners,'planecorners',planecorners);
    geofiledata.firstcornertoskip = 4;
    Sindata = struct('coordinates',[-1 -1 0]);
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
    Rindata = struct('coordinates',receivers);
    controlparameters = struct('frequencies',100);
    controlparameters.difforder = 4;
    controlparameters.ngauss = 32;
    filehandlingparameters = struct('outputdirectory',infilepath);
    filehandlingparameters.filestem = filestem;
    filehandlingparameters.savelogfile = 1;
    filehandlingparameters.showtext = 1;
    controlparameters.Rstart = 0;

    envdata.cair = 344;

    EDmain_convexESIE(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);

    eval(['load ',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tf.mat tfdirect tfgeom tfdiff EDversionnumber'])
    eval(['load ',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tfinteq.mat'])

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

    fwrite(fid,[' ',lineending],'char');
    fwrite(fid,['####################################################################',lineending],'char');
    fwrite(fid,['Test ',II,': EDmain_convexESIE, diff1 continuity across corner zone boundary, at 100 Hz',lineending],'char');
    fwrite(fid,['Two edges meat at 90 deg.; three receivers exactly at, and very near zone boundaries',lineending],'char');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Test 4: EDmain_convexESIE, diff1 fine-grain continuity across zone boundary, at 100 Hz
%

if runtest(4) == 1
    itest = 4;
    II = int2str(itest);
    
    if showtext > 0
        disp(' ')
        disp('*********************************************************************')
        disp(['Test ',II,': EDmain_convexESIE, diff1 fine-grain continuity across zone boundary, at 100 Hz']);
        disp('Single edge; three receivers exactly at, and very near zone boundaries')
        disp('Computed value at ZB should be very close to the mean value of the two')
        disp('surrounding responses (< 1e-5)')
    else
        disp(['Test ',II,': EDmain_convexESIE, diff1 continuity across zone boundary, at 100 Hz']);        
    end
    mfile = mfilename('fullpath');
    [infilepath,filestem] = fileparts(mfile);

    corners = [     0   0   1
        0 0 -1
        -1 0 -1
        -1 0 1];

    planecorners = [4 3 2 1 ;1 2 3 4];

    geofiledata = struct('corners',corners,'planecorners',planecorners);
    % geofiledata.firstcornertoskip = 3;
    Sindata = struct('coordinates',[-1 1 0]/sqrt(2));
    phivecd = [45.05 45.01 45.001 45.0 44.999 44.99 44.95 315.05 315.01 315.001 315.0 314.999 314.99 314.95].';
    receivers = [cosd(phivecd) sind(phivecd) 0*phivecd];
    Rindata = struct('coordinates',receivers);
    controlparameters = struct('frequencies',100);
    controlparameters.difforder = 1;
    filehandlingparameters = struct('outputdirectory',infilepath);
    filehandlingparameters.filestem = filestem;
    filehandlingparameters.savelogfile = 1;
    filehandlingparameters.savelogfile = 1;
    filehandlingparameters.showtext = 0;
    controlparameters.Rstart = 0;

    envdata.cair = 344;

    EDmain_convexESIE(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);

    eval(['load ',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tf.mat tfdirect tfgeom tfdiff EDversionnumber'])

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
       title(['Test ',II,' Direct sound amplitude across zone boundary']);
       figure
       plot(phivecd(iv2),abs(tfZBspec),'-o',phivecd(iv2(ivmid2)),abs(tfZBspec(ivmid2)),'*')
       xlabel('Receiver angle   [deg.]');
       ylabel('TF magnitude   [-]');
       title(['Test ',II,' Specular reflection amplitude across zone boundary']);
       
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fwrite(fid,[' ',lineending],'char');
    fwrite(fid,['####################################################################',lineending],'char');
    fwrite(fid,['Test ',II,': EDmain_convexESIE, diff1 fine-grain continuity across zone boundary, at 100 Hz',lineending],'char');
    fwrite(fid,['Single edge; seven receivers exactly at, and very near zone boundaries',lineending],'char');
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose(fid);






