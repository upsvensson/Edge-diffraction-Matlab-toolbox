% EDexample_sqplateD10_circpistonnearfield.m
%

mfile = mfilename('fullpath'); 
[infilepath,filestem] = fileparts(mfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A 60cm*60cm dummy plate is generated, acting as an infinite baffle

corners = 0.3*[1 1 0;-1 1 0;-1 -1 0; 1 -1 0];
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

fvec = [50 ].';
controlparameters = struct('frequencies',fvec); 
controlparameters.docalctf = 1;
controlparameters.difforder = 10;
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
tfGA_piston = EDres.tfdirect;
tfD1_piston = EDres.tfdirect + EDres.tfdiff;

figure(2)
clf(2)
h = semilogx(zvec,20*log10(abs([tfGA_piston ])),'o',zvec,20*log10(abs([tfD1_piston ])),'o',...
    zvec,20*log10(abs([tftot_piston ])),'*',zvec,20*log10(2./zvec),'-',zvec,20*log10(1./zvec),'-'); 
g = get(h(1),'Parent');
set(g,'FontSize',14);
ncurves = length(h);
for ii = 1:ncurves
    set(h(ii),'LineWidth',2);
end
g = xlabel('On-axis distance [m]');
set(g,'FontSize',14)
g = ylabel('TF magnitude re. on-axis [dB]');
set(g,'FontSize',14)
g = title('On-axis response for a 20 cm piston at 100 Hz'); 
set(g,'FontSize',14)
%xlim([50 5000])
grid
g = legend('Piston GA (Rayleigh)','Piston w D1','Piston w D10','2/distance','1/distance');
set(g,'Location','best')
set(g,'FontSize',14)

figure(3)
clf(3)
h = semilogx(zvec,20*log10(abs([tfD1_piston ]))-20*log10(1./zvec),'-o',...
    zvec,20*log10(abs([tftot_piston ]))-20*log10(1./zvec),'-*');
g = get(h(1),'Parent');
set(g,'FontSize',14);
ncurves = length(h);
for ii = 1:ncurves
    set(h(ii),'LineWidth',2);
end
g = xlabel('On-axis distance [m]');
set(g,'FontSize',14)
g = ylabel('TF magnitude re. on-axis [dB]');
set(g,'FontSize',14)
g = title('On-axis response for a 20 cm piston at 50 Hz'); 
set(g,'FontSize',14)
%xlim([50 5000])
grid
g = legend('Piston w D1 re 1/distance','Piston w D10 re 1/distance');
set(g,'Location','best')
set(g,'FontSize',14)

figure(1)
clf(1)
eddatafile = [infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_eddata.mat']; 
EDplotmodel(eddatafile,3);
