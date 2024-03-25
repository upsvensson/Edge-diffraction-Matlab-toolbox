% EDexample_LspKessel_minimal.m

mfile = mfilename('fullpath'); 
[infilepath,filestem] = fileparts(mfile);

%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
% Here the size of the shoebox−shaped enclosure is defined
% The function EDmakegeo_shoebox centers the box around the origin in the
% xy-direction but one face is at z = 0 and the opposite face at -l_z.
% We shift it in the y−direction so that the origin (where the source will
% be) has 0.2m to all edges.

[corners,planecorners] = EDmakegeo_shoebox(0.4,0.64,0.32); 
corners(:,2) = corners(:,2) - 0.12;
geoinputdata = struct('corners',corners,'planecorners',planecorners);

%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−− 
% The source (single point) and receiver are defined

sourcecoordinates = [0 0 0.00001];
Sinputdata = struct('coordinates',sourcecoordinates);

receivercoordinates = [0 0 1];
Rinputdata = struct('coordinates',receivercoordinates);

%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−− 
% Set some calculation parameter values

fvec = linspace(50,3000,100);
controlparameters = struct('frequencies',fvec); 
controlparameters.docalctf = 1;

% The settings below could be changed, if desired
% .ngauss could be increased for higher frequencies, and/or higher accuracy 
% .difforder could be reduced for faster calculations (the default is 10)
% controlparameters.ngauss = 24;
% controlparameters.difforder = 6;

%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−− 
% Specify output file names and location
% If the directory doesn't exist; it will be created

outputdirectory = [infilepath,filesep,'results'];
filehandlingparameters = struct('outputdirectory',outputdirectory); 
filehandlingparameters.filestem = filestem;

%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−− 
% Run the calculations

EDres = EDmain_convex(geoinputdata,Sinputdata,... 
   Rinputdata,struct,controlparameters,filehandlingparameters);

%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
% Present the results

figure(1) 
semilogx(controlparameters.frequencies,20*log10(abs(EDres.tftot)),'-o') 
xlabel('Frequency [Hz]');
ylabel('TF magnitude re. 1m free−field [dB]')
title('Frequency response of the Kessel loudspeaker, at 1m distance') 
axis([50 5000 0 10]);
grid

EDplotmodel([filehandlingparameters.outputdirectory,filesep,... 
    filehandlingparameters.filestem,'_eddata.mat'],'plotoptions',3,'figurewindow',2);
