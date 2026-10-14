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

Sdist = 1e6;
sourcecoordinates = Sdist/sqrt(3)*[1 1 1];
Sinputdata = struct('coordinates',sourcecoordinates);

receivercoordinates = [0 0 0.00001];
Rinputdata = struct('coordinates',receivercoordinates);

%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−− 
% Set some calculation parameter values

fvec = 0; %linspace(50,3000,100);
controlparameters = struct('frequencies',fvec); 
controlparameters.docalctf = 1;

% The settings below could be changed, if desired
% .ngauss could be increased for higher frequencies, and/or higher accuracy 
% .difforder could be reduced for faster calculations (the default is 10)
%controlparameters.ngauss = 24;
%ngaussvec = [24 36 48 60 72 84 96 124 148];
ngaussvec = [24 36 48 60 72 84 96];
ncases = length(ngaussvec);
controlparameters.difforder = 15;

%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−− 
% Specify output file names and location
% If the directory doesn't exist; it will be created

outputdirectory = [infilepath,filesep,'results'];
filehandlingparameters = struct('outputdirectory',outputdirectory); 
filehandlingparameters.filestem = filestem;

%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−− 
% Run the calculations

allres = zeros(1,ncases);
for ii = 1:ncases
    controlparameters.ngauss = ngaussvec(ii);
    filehandlingparameters.filestem = [filestem,int2str(ii)];
    EDres = EDmain_convex(geoinputdata,Sinputdata,... 
    Rinputdata,struct,controlparameters,filehandlingparameters);
    allres(ii) = EDres.tftot;
end


extrapolresult = 0*allres;
allxvalues = 1./ngaussvec.^2;
allyvalues = allres*Sdist;
for ii = 3:ncases
    p = polyfit(allxvalues(ii-2:ii),allyvalues(ii-2:ii),1);
    extrapolresult(ii) = p(2);
end

figure(1)
clf(1)
h = plot(1./ngaussvec.^2,allres*Sdist,'-o',allxvalues(3:end),extrapolresult(3:end),'-*');
seth
sethlines
grid



return

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
