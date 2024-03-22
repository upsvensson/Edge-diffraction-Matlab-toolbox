% EDexample_Rayleigh_circpistondirectivity.m
%

mfile = mfilename('fullpath'); 
[infilepath,filestem] = fileparts(mfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A 10m*10m dummy plate is generated, acting as an infinite baffle

corners = 1*[1 1 0;-1 1 0;-1 -1 0; 1 -1 0];
planecorners = [1 2 3 4;4 3 2 1];
geoinputdata = struct('corners',corners,'planecorners',planecorners);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The source is modeled as a 32-sided regular polygon approximation of a
% circular piston. The "polygonradius" is adjusted so that the regular
% polygon gets the same area as the circular piston that it tries to
% represent.
%
% To study the directivity, we create an arc of receiver

circpistonradius = 0.1;
npistoncorners = 32;
polygonradius = circpistonradius*sqrt(2*pi/npistoncorners/sin(2*pi/npistoncorners));
phivec = [0:npistoncorners-1].'*2*pi/npistoncorners;
pistoncorners = polygonradius*[cos(phivec) sin(phivec) zeros(size(phivec))];

Sinputdata = struct('sourcetype','polygonpiston');
Sinputdata.pistoncornercoordinates = pistoncorners;
Sinputdata.pistoncornernumbers = [1:npistoncorners];
Sinputdata.pistonplanes = 1;

thetavec = linspace(0,0.999*pi/2,100);

Rdistance = 100;
receivercoordinates = Rdistance*[sin(thetavec.') zeros(size(thetavec.')) cos(thetavec.')] ;
Rinputdata = struct('coordinates',receivercoordinates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Some calculation parameters
% No diffraction.
% We choose 3 frequencies to illustrate three degrees of directivity

fvec = [100 1e3 5e3].';
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

tftot_piston = EDres.tftot*Rdistance/2;
cair = 343;
kavec = 2*pi*fvec/cair*circpistonradius;

figure(2)
clf(2)
h = plot(thetavec*180/pi,(abs([tftot_piston ])),'-o'); 
g = get(h(1),'Parent');
set(g,'FontSize',14);
set(h(1),'LineWidth',2);
g = xlabel('Radiation angle [deg.]');
set(g,'FontSize',14)
g = ylabel('TF magnitude re. on-axis [-]');
set(g,'FontSize',14)
g = title('Directivity function of a 20 cm piston'); 
set(g,'FontSize',14)
%xlim([50 5000])
grid
g = legend(['100 Hz (ka=',num2str(kavec(1)),')'],...
           ['1 kHz (ka=',num2str(kavec(2)),')'],...
           ['5 kHz (ka=',num2str(kavec(3)),')']);
set(g,'Location','best')
set(g,'FontSize',14)

figure(1)
clf(1)
eddatafile = [infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_eddata.mat']; 
EDplotmodel(eddatafile,'plotoption',3,'figurewindow',1);

figure(3)
clf(3)
EDplotmodel(eddatafile,'plotoption',1,'figurewindow',3);
