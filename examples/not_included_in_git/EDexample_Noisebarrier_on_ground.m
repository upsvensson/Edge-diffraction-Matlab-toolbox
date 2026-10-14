% EDexample_Noisebarrier_no_ground.m
%

mfile = mfilename('fullpath'); 
[infilepath,filestem] = fileparts(mfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The noise barrier is created as a shoebox, with double the height of the
% physical barrier, to form an image copy which represents the rigid
% ground.
% The function EDmakegeo_shoebox centers the box around the origin in the
% xy-direction but one face is at z = 0 and the opposite face at -l_z.

noisebarrierheight = 3;
noisebarrierthickness = 0.5;
noisebarrierlength = 50;
[corners,planecorners] = EDmakegeo_shoebox(noisebarrierthickness,noisebarrierlength,2*noisebarrierheight); 
corners(:,3) = corners(:,3) + noisebarrierheight;
geoinputdata = struct('corners',corners,'planecorners',planecorners);
geoinputdata.listofedgestoskip = [1 2;3 4;5 6;7 8];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The source will be two monopoles (the physical source and its image
% version), and their contributions should be added.

sourcecoordinates = [-3 0 1];
sourcecoordinates = [sourcecoordinates;sourcecoordinates];
sourcecoordinates(2,3) = -sourcecoordinates(2,3);
Sinputdata = struct('coordinates',sourcecoordinates);
Sinputdata.doaddsources = 1;

receivercoordinates = [2 0 1];
Rinputdata = struct('coordinates',receivercoordinates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Some calculation parameters

fvec = linspace(50,500,300);
controlparameters = struct('frequencies',fvec); 
controlparameters.docalctf = 1;
controlparameters.difforder = 4;
controlparameters.ngauss = 32;
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

tftot = EDres.tftot; 

figure(2)
clf(2)
h = semilogx(controlparameters.frequencies,20*log10(abs(tftot)),'-o'); 
g = get(h(1),'Parent');
set(g,'FontSize',14);
set(h(1),'LineWidth',2);
g = xlabel('Frequency [Hz]');
set(g,'FontSize',14)
g = ylabel('TF magnitude re. 1m free-field [dB]');
set(g,'FontSize',14)
g = title('Freq. resp. of a noise barrier on a rigid ground'); 
set(g,'FontSize',14)
xlim([50 5000])
grid

figure(1)
clf(1)
eddatafile = [infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_eddata.mat']; 
EDplotmodel(eddatafile,3,'figurewindow',1);

figure(3)
clf(3)
eddatafile = [infilepath,filesep,'results',filesep,filehandlingparameters_monopole.filestem,'_eddata.mat']; 
EDplotmodel(eddatafile,3,'figurewindow',3,'offedges',1);
