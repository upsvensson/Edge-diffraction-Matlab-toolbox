% EDexample_LspCylinder16.m 
%

mfile = mfilename('fullpath'); 
[infilepath,filestem] = fileparts(mfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the cylinder geometry using the function EDmakegeo_cylinder

radius = 0.3048;
length = 2*radius;
[corners,planecorners,ncorners,radius] = EDmakegeo_cylinder(radius,length,16,'e',1,0,-radius); 
geoinputdata = struct('corners',corners,'planecorners',planecorners);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% The source (single point) and receiver are defined

sourcecoordinates = [0 0 0.00001];
Sinputdata = struct('coordinates',sourcecoordinates);

receivercoordinates = [0 0 6.6];
Rinputdata = struct('coordinates',receivercoordinates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Some calculation parameters

fvec = linspace(50,4000,100);
controlparameters = struct('frequencies',fvec); 
controlparameters.docalctf = 1; 
controlparameters.ngauss = 24; 
controlparameters.difforder = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Output file names and location

filehandlingparameters = struct('outputdirectory',[infilepath,filesep,'results']); 
filehandlingparameters.filestem = filestem;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Run the calculations

EDres = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Present the results

figure(2)
clf(2)
h = semilogx(controlparameters.frequencies,20*log10(abs(EDres.tftot)),'−o');
g = get(h(1),'Parent');
set(g,'FontSize',14);
set(h(1),'LineWidth',2);
g = xlabel('Frequency [Hz]');
set(g,'FontSize',14)
g = ylabel('TF magnitude re. 1m   [dB]')
set(g,'FontSize',14)
g = title('Frequency response of the Kessel loudspeaker, at 1m distance, with a 20 cm piston'); 
set(g,'FontSize',14)
xlim([50 5000])
grid

figure(1)
clf(1)
eddatafile = [infilepath,filesep,'results',filesep,filehandlingparameters_piston.filestem,'_eddata.mat']; 
EDplotmodel(eddatafile,1)