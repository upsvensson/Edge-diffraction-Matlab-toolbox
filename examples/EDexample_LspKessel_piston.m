% EDexample_LspKessel_piston.m
%

mfile = mfilename('fullpath'); 
[infilepath,filestem] = fileparts(mfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here the size of the shoebox−shaped enclosure is defined 
% The function EDmakegeo_shoebox centers the box around the origin in the
% xy-direction but one face is at z = 0 and the opposite face at -l_z.
% We shift it in the y−direction so that the origin (where the source will
% be) has 0.2m to all edges.

[corners,planecorners] = EDmakegeo_shoebox(0.4,0.64,0.32); 
corners(:,2) = corners(:,2) - 0.12;
geoinputdata = struct('corners',corners,'planecorners',planecorners);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The source (piston−like distribution of monopoles, using the function EDcirclepoints) and receiver are defined

sourcecoordinates = EDcirclepoints(0.1,2); sourcecoordinates(:,3) = 0.00001;
nsources = size(sourcecoordinates,1);
Sinputdata = struct('coordinates',sourcecoordinates); 
Sinputdata.doaddsources = 1;
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

filehandlingparameters_piston.filestem = filehandlingparameters.filestem;
EDres_piston = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the changes for a single monopole source and rerun the calculations

Sinputdata.coordinates = Sinputdata.coordinates(1,:);
filehandlingparameters.filestem = [filehandlingparameters.filestem,'_monopole'];
EDres_monopole = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Present the results

tftot_piston = EDres_piston.tftot/nsources; 
tftot_point = EDres_monopole.tftot;

figure(2)
clf(2)
h = semilogx(controlparameters.frequencies,20*log10(abs([tftot_piston tftot_point])),'-o'); 
g = get(h(1),'Parent');
set(g,'FontSize',14);
set(h(1),'LineWidth',2);
set(h(2),'LineWidth',2);
g = xlabel('Frequency [Hz]');
set(g,'FontSize',14)
g = ylabel('TF magnitude re. 1m [dB]');
set(g,'FontSize',14)
g = title('Frequency response of the Kessel loudspeaker, at 1m distance, with a 20 cm piston'); 
set(g,'FontSize',14)
xlim([50 5000])
grid
g = legend('Piston source','Point source');
set(g,'Location','best')
set(g,'FontSize',14)

figure(1)
clf(1)
eddatafile = [infilepath,filesep,'results',filesep,filehandlingparameters_piston.filestem,'_eddata.mat']; 
EDplotmodel(eddatafile,1)
