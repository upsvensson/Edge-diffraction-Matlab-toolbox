% EDexample_Rayleigh_circpistonnearfield.m
%

mfile = mfilename('fullpath'); 
[infilepath,filestem] = fileparts(mfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A 2m*2m dummy plate is generated, acting as an infinite baffle

corners = 1*[1 1 0;-1 1 0;-1 -1 0; 1 -1 0];
planecorners = [1 2 3 4;4 3 2 1];
geoinputdata = struct('corners',corners,'planecorners',planecorners);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The source is modeled as a 32-sided regular polygon approximation of a
% circular piston. The "polygonradius" is adjusted so that the regular
% polygon gets the same area as the circular piston that it tries to
% represent.
%
% The receivers are placed along the on-axis

circpistonradius = 0.1;
npistoncorners = 32;
polygonradius = circpistonradius*sqrt(2*pi/npistoncorners/sin(2*pi/npistoncorners));
phivec = [0:npistoncorners-1].'*2*pi/npistoncorners;
pistoncorners = polygonradius*[cos(phivec) sin(phivec) zeros(size(phivec))];

Sinputdata = struct('sourcetype','polygonpiston');
Sinputdata.pistoncornercoordinates = pistoncorners;
Sinputdata.pistoncornernumbers = [1:npistoncorners];
Sinputdata.pistonplanes = 1;

nreceivers = 100;
zvec = logspace(log10(1e-3),log10(10),nreceivers);
receivercoordinates = [zeros(nreceivers,2) zvec(:)];
Rinputdata = struct('coordinates',receivercoordinates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Some calculation parameters
% No diffraction.
% We choose 1 frequency

fvec = [100 ].';
controlparameters = struct('frequencies',fvec); 
controlparameters.docalctf = 1;
controlparameters.difforder = 0;
nfrequencies = length(controlparameters.frequencies);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Output file names and location

filehandlingparameters = struct('outputdirectory',[infilepath,filesep,'results']); 
filehandlingparameters.filestem = filestem;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Run the calculations

EDres = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Present the results

tftot_piston = EDres.tftot;

figure(2)
clf(2)
h = loglog(zvec,(abs([tftot_piston ])),'-o',zvec,2./zvec,'-'); 
g = get(h(1),'Parent');
set(g,'FontSize',14);
set(h(1),'LineWidth',2);
set(h(2),'LineWidth',2);
g = xlabel('On-axis distance [m]');
set(g,'FontSize',14)
g = ylabel('TF magnitude re. on-axis [-]');
set(g,'FontSize',14)
g = title('On-axis response for a 20 cm piston at 100 Hz'); 
set(g,'FontSize',14)
%xlim([50 5000])
grid
g = legend('Piston response','2/distance');
set(g,'Location','best')
set(g,'FontSize',14)

figure(1)
clf(1)
eddatafile = [infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_eddata.mat']; 
EDplotmodel(eddatafile,3);
