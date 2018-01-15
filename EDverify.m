% EDverify.m
%
% Script for verifying the EDtoolbox.
% It contains a number of tests which test various aspects of the code
% 
% Peter Svensson 15 Jan. 2018 (peter.svensson@ntnu.no)

% 15 Jan. 2018 First version

showtext_verify = 0;

ntests = 3;
passtest = zeros(ntests,1);

[EDversionnumber,changedate,changetime] = EDversion;

clockvec = clock;
datevec = date;
datevec = datevec(datevec~='-');
datetimevec = [datevec,'_',sprintf('%02d',clockvec(4)),'h', sprintf('%02d',clockvec(5)),'m',sprintf('%02d',round(clockvec(6)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Test 1
%
% Response at surface, plane wave incidence, DC, no singularity

if showtext_verify > 0
    disp(' ')
    disp(['EDverify, EDtoolbox v. ',num2str(EDversionnumber),' (last change on ',changedate,'), run on ',datetimevec])
    disp(' ')
    disp('*********************************************************************')
    disp('Test 1: EDmain_convexESIE, DC response at surface, plane wave incidence');
    disp('Expected rel. errors for receiver positions 1 (in front, at surface) and 2 (behind, at surface), with')
    disp('Gauss-Legendre quadrature, and 96 points for the longest edge, and diffraction order 20, are')
    disp('1.3e-5 and 1.5e-5, respectively.')
else
    disp(' ')
    disp(['EDverify, EDtoolbox v. ',num2str(EDversionnumber),' (last change on ',changedate,'), run on ',datetimevec])
    disp(' ')
    disp('*********************************************************************')
    disp('Test 1: EDmain_convexESIE, DC response at surface, plane wave incidence');    
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
eval(['load ',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tf.mat tfdirect tfgeom tfdiff EDversion'])
    
tftot = tfdirect + tfgeom + tfdiff + tfinteqdiff;
tftot = tftot*norm(Sindata.coordinates);

relerr = abs(tftot-1);
if relerr(1) < 0.14e-4 && relerr(2) < 0.16e-4
    passtest(1) = 1;
end
    
if showtext_verify > 0
    disp(' ')
    disp(['Computed results are:   ',num2str(relerr(1)),' and ',num2str(relerr(2))])
    disp(' ')
    if passtest(1) == 1
        disp('So, verification test 1 was passed')
    else
        disp('So, verification test 1 was not passed. Please check the code')
        disp('Make sure these settings were set in this script: difforder = 2, ngauss = 96, discretizationtype = 2')
    end
end
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

logfilename = [infilepath,filesep,'results',filesep,'EDverify_v',num2str(EDversion),'_',datetimevec,'.txt'];

fid = fopen(logfilename,'w');
if fid == -1
    disp('The logfile is not possible to open - check that it isn''t opened by any program!')
    return
end
fwrite(fid,['EDverify, EDtoolbox v. ',num2str(EDversionnumber),' (last change on ',changedate,'), run on ',datetimevec,lineending],'char');
fwrite(fid,[' ',lineending],'char');
fwrite(fid,['####################################################################',lineending],'char');
fwrite(fid,['Test 1: EDmain_convexESIE, DC response at surface, plane wave incidence',lineending],'char');
fwrite(fid,['Expected rel. errors for receiver positions 1 (in front, at surface) and 2 (behind, at surface), with',lineending],'char');
fwrite(fid,['Gauss-Legendre quadrature, and 96 points for the longest edge, and diffraction order 20, are',lineending],'char');
fwrite(fid,['1.3e-5 and 1.5e-5, respectively.',lineending],'char');
fwrite(fid,[' ',lineending],'char');
fwrite(fid,['Computed results are: ',num2str(relerr(1)),' and ',num2str(relerr(2)),lineending],'char');
if relerr(1) < 0.14e-4 && relerr(2) < 0.16e-4
    fwrite(fid,['So, verification test 1 was passed',lineending],'char');
else
    fwrite(fid,['So, verification test 1 was not passed. Please check the code'   ,lineending],'char');
    fwrite(fid,['Make sure these settings were set in this script: difforder = 2, ngauss = 96, discretizationtype = 2',lineending],'char');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Test 2: EDmain_convexESIE, diff1 continuity across zone boundary, at 100 Hz
%

if showtext_verify > 0
    disp(' ')
    disp('*********************************************************************')
    disp('Test 2: EDmain_convexESIE, diff1 continuity across zone boundary, at 100 Hz');
    disp('Single edge; three receivers exactly at, and very near zone boundaries')
    disp('Computed value at ZB should be very close to the mean value of the two')
    disp('surrounding responses (< 1e-5)')
else
    disp('Test 2: EDmain_convexESIE, diff1 continuity across zone boundary, at 100 Hz');        
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
    
eval(['load ',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tf.mat tfdirect tfgeom tfdiff EDversion'])
    
tftot = tfdirect + tfgeom + tfdiff;

tfZBdirect = tfdirect(4:6) + tfdiff(4:6);
tfZBdirectmean = (tfZBdirect(1) + tfZBdirect(3))/2;
relerrZBdirect = abs( (tfZBdirectmean-tfZBdirect(2))/tfZBdirectmean );

tfZBspec = tfgeom(1:3) + tfdiff(1:3);
tfZBspecmean = (tfZBspec(1) + tfZBspec(3))/2;
relerrZBspec = abs( (tfZBspecmean-tfZBspec(2))/tfZBspecmean );

if relerrZBdirect < 1e-5 & relerrZBspec < 1e-5
   passtest(2) = 1; 
end

if showtext_verify > 0
    disp(' ')
    disp(['Computed results at the direct sound ZB differ by: ',num2str(relerrZBdirect)])
    disp(['Computed results at the specular reflection ZB differ by: ',num2str(relerrZBspec)])
    disp(' ')
    if passtest(2) == 1
        disp('So, verification test 2 was passed')
    else
        disp('So, verification test 2 was not passed. Please check the code')   
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fwrite(fid,[' ',lineending],'char');
fwrite(fid,['####################################################################',lineending],'char');
fwrite(fid,['Test 2: EDmain_convexESIE, diff1 continuity across zone boundary, at 100 Hz',lineending],'char');
fwrite(fid,['Single edge; three receivers exactly at, and very near zone boundaries',lineending],'char');
fwrite(fid,['Computed value at ZB should be very close to the mean value of the two',lineending],'char');
fwrite(fid,['surrounding responses (< 1e-5)',lineending],'char');
fwrite(fid,[' ',lineending],'char');
fwrite(fid,['Computed results at the direct sound ZB differ by: ',num2str(relerrZBdirect),lineending],'char');
fwrite(fid,['Computed results at the specular reflection ZB differ by: ',num2str(relerrZBspec),lineending],'char');
fwrite(fid,[' ',lineending],'char');
if relerrZBdirect < 1e-5 & relerrZBspec < 1e-5
    fwrite(fid,['So, verification test 2 was passed',lineending],'char');
else
    fwrite(fid,['So, verification test 2 was not passed. Please check the code'   ,lineending],'char');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Test 3: EDmain_convexESIE, diff1 continuity across corner zone boundary, at 100 Hz
%

if showtext_verify > 0
    disp(' ')
    disp('*********************************************************************')
    disp('Test 3: EDmain_convexESIE, diff1 continuity across corner zone boundary, at 100 Hz');
    disp('Two edges meat at 90 deg.; three receivers exactly at, and very near zone boundaries')
    disp('Computed value at ZB should be very close to the mean value of the two')
    disp('surrounding responses (< 1e-3)')
else
    disp('Test 3: EDmain_convexESIE, diff1 continuity across corner zone boundary, at 100 Hz');    
end

mfile = mfilename('fullpath');
[infilepath,filestem] = fileparts(mfile);

corners = [     0   0   -1
    0 0 1
    -1 0 1
    -1 0 -1];

planecorners = [4 3 2 1 ;1 2 3 4];

geofiledata = struct('corners',corners,'planecorners',planecorners);
geofiledata.firstcornertoskip = 4;
Sindata = struct('coordinates',[-1 -1 0]);
receivers =[0.995 -1 1.99; 
            1 -1 2; 
            1.005 -1 2.01;
            0.995 1 1.99;
            1 1 2; 
            1.005 1 2.01];
Rindata = struct('coordinates',receivers);
controlparameters = struct('frequencies',100);
controlparameters.difforder = 4;
controlparameters.ngauss = 32;
filehandlingparameters = struct('outputdirectory',infilepath);
filehandlingparameters.filestem = filestem;
filehandlingparameters.savelogfile = 1;
filehandlingparameters.showtext = 0;
controlparameters.Rstart = 0;

envdata.cair = 344;

EDmain_convexESIE(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);
    
eval(['load ',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tf.mat tfdirect tfgeom tfdiff EDversion'])
eval(['load ',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tfinteq.mat'])
    
tftot = tfdirect + tfgeom + tfdiff;

tfZBdirect = tfdirect(4:6) + tfdiff(4:6) + tfinteqdiff(4:6);
tfZBdirectmean = (tfZBdirect(1) + tfZBdirect(3))/2;
relerrZBdirect = abs( (tfZBdirectmean-tfZBdirect(2))/tfZBdirectmean );

tfZBspec = tfgeom(1:3) + tfdiff(1:3) + tfinteqdiff(1:3);
tfZBspecmean = (tfZBspec(1) + tfZBspec(3))/2;
relerrZBspec = abs( (tfZBspecmean-tfZBspec(2))/tfZBspecmean );

if relerrZBdirect < 1e-3 & relerrZBspec < 1e-3
   passtest(3) = 1; 
end

if showtext_verify > 0
     disp(' ')
    disp(['Computed results at the direct sound corner ZB differ by: ',num2str(relerrZBdirect)])
    disp(['Computed results at the specular reflection corner ZB differ by: ',num2str(relerrZBspec)])
    disp(['   '])
    if passtest(3) == 1
        disp('So, verification test 3 was passed')
    else
        disp('So, verification test 3 was not passed. Please check the code')   
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fwrite(fid,[' ',lineending],'char');
fwrite(fid,['####################################################################',lineending],'char');
fwrite(fid,['Test 3: EDmain_convexESIE, diff1 continuity across corner zone boundary, at 100 Hz',lineending],'char');
fwrite(fid,['Two edges meat at 90 deg.; three receivers exactly at, and very near zone boundaries',lineending],'char');
fwrite(fid,['Computed value at ZB should be very close to the mean value of the two',lineending],'char');
fwrite(fid,['surrounding responses (< 1e-3)',lineending],'char');
fwrite(fid,[' ',lineending],'char');
fwrite(fid,['Computed results at the direct sound corner ZB differ by: ',num2str(relerrZBdirect),lineending],'char');
fwrite(fid,['Computed results at the specular reflection corner ZB differ by: ',num2str(relerrZBspec),lineending],'char');
fwrite(fid,[' ',lineending],'char');
if relerrZBdirect < 1e-3 & relerrZBspec < 1e-3
    fwrite(fid,['So, verification test 3 was passed',lineending],'char');
else
    fwrite(fid,['So, verification test 3 was not passed. Please check the code'   ,lineending],'char');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose(fid);






