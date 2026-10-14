% EDexample_Cube_surfacepressure.m

mfile = mfilename('fullpath');
[infilepath,filestem] = fileparts(mfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the cube
[corners,planecorners] = EDmakegeo_shoebox(1,1,1);
corners(:,3) = corners(:,3) + 0.5;
geoinputdata = struct('corners',corners,'planecorners',planecorners);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the source coordinate

Sdist = 1e5;
Scoordinates = Sdist*[1 1 1]/sqrt(3);
Sinputdata = struct('coordinates',Scoordinates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the receiver coordinate

xvec = linspace(-0.49,0.49,10);
yvec = linspace(-0.49,0.49,10);

Rcoordinates = [xvec(:),0*xvec(:),0*xvec(:) + 0.5001];

Rinputdata = struct('coordinates',Rcoordinates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation settings

controlparameters = struct('difforder',10);
controlparameters.frequencies = [10 20 40 80 160];
controlparameters.docalctf = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify where result files will be saved

filehandlingparameters = struct('outputdirectory',[infilepath,filesep,'results']);
filehandlingparameters.filestem = filestem;

fvec = controlparameters.frequencies;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation settings, for the TF calculations

controlparameters.ngauss = 24;
controlparameters.difforder = 10;
controlparameters.docalctf = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the TF calculations

EDres = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Present the results
    
tftot = EDres.tftot*Sdist;

figure(2)
clf(2)
h = semilogx(fvec,20*log10(abs(tftot)),'*');
for ii = 1:4
   set(h(ii),'LineWidth',2) 
end
g = get(h(1),'Parent');
set(g,'FontSize',14)
grid
g = xlabel('Frequency   [Hz]');
set(g,'FontSize',14)
g = ylabel('TF magnitude re. 1m   [dB]');
set(g,'FontSize',14)
g = legend('TF from ir: GA only','TF from ir incl. diffr.1','TF from ir incl. diffr.2','TF from ir incl. diffr.3',['FD TF: Incl. diffr.',int2str(controlparameters.difforder)]);
set(g,'FontSize',14,'Location','SouthEast')
%axis([20 20000 -8 4])

figure(3)
clf(3)
eddatafile = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_eddata.mat'];
EDplotmodel(eddatafile,'plotoptions',3,'figurewindow',3);

