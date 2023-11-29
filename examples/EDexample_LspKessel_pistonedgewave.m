% EDexample_LspKessel_pistonedgewave.m
%

mfile = mfilename('fullpath'); 
[infilepath,filestem] = fileparts(mfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here the size of the shoebox−shaped enclosure is defined 
% The function EDmakegeo_shoebox centers the box around the origin but we
% shift it in the x−direction so that the origin (where the source will be
% has 0.2m to all edges)

[corners,planecorners] = EDmakegeo_shoebox(0.4,0.64,0.32); 
corners(:,2) = corners(:,2) - 0.12;
geoinputdata = struct('corners',corners,'planecorners',planecorners);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The source is modeled as a 32-sided regular polygon approximation of a
% circular piston. The "polygonradius" is adjusted so that the regular
% polygon gets the same area as the circular piston that it tries to
% represent.
% For the diffraction terms, a number of monopoles are representing the
% piston.

circpistonradius = 0.1;
npistoncorners = 32;
polygonradius = circpistonradius*sqrt(2*pi/npistoncorners/sin(2*pi/npistoncorners));
phivec = [0:npistoncorners-1].'*2*pi/npistoncorners;
pistoncorners = polygonradius*[cos(phivec) sin(phivec) zeros(size(phivec))];

ngaussorder = 3;

Sinputdata = struct('sourcetype','polygonpiston');
Sinputdata.pistoncornercoordinates = pistoncorners;
Sinputdata.pistoncornernumbers = [1:npistoncorners];
Sinputdata.pistonplanes = 1;
Sinputdata.pistongaussorder = ngaussorder;

receivercoordinates = [0 0 1];
Rinputdata = struct('coordinates',receivercoordinates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Some calculation parameters

fvec = linspace(50,3000,100);
controlparameters = struct('frequencies',fvec); 
controlparameters.docalctf = 1;
nfrequencies = length(controlparameters.frequencies);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Output file names and location

filehandlingparameters = struct('outputdirectory',[infilepath,filesep,'results']); 
filehandlingparameters.filestem = filestem;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Run the calculations, for the piston case

EDres = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Present the results

tftot_piston = EDres.tftot; 

figure(2)
clf(2)
h = semilogx(controlparameters.frequencies,20*log10(abs([tftot_piston ])),'-o'); 
g = get(h(1),'Parent');
set(g,'FontSize',14);
set(h(1),'LineWidth',2);
g = xlabel('Frequency [Hz]');
set(g,'FontSize',14)
g = ylabel('TF magnitude re. 1m [dB]');
set(g,'FontSize',14)
g = title('Frequency response of the Kessel loudspeaker, at 1m distance, with a 20 cm piston'); 
set(g,'FontSize',14)
xlim([50 5000])
grid
g = legend('Piston source');
set(g,'Location','best')
set(g,'FontSize',14)

figure(1)
clf(1)
eddatafile = [infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_eddata.mat']; 
EDplotmodel(eddatafile,3);
