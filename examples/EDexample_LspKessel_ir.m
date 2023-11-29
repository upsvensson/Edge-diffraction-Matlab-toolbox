% EDexample_LspKessel_ir.m

mfile = mfilename('fullpath');
[infilepath,filestem] = fileparts(mfile);

[corners,planecorners] = EDmakegeo_shoebox(0.4,0.64,0.32);
corners(:,2) = corners(:,2) - 0.12;
geoinputdata = struct('corners',corners,'planecorners',planecorners);
Sinputdata = struct('coordinates',[0 0 0.0001]);
Rinputdata = struct('coordinates',[0 0 2]);
controlparameters = struct('fs',48000);
controlparameters.docalcir = 1;
controlparameters.difforder = 3;
controlparameters.savealldifforders = 1;
filehandlingparameters = struct('outputdirectory',[infilepath,filesep,'results']);
filehandlingparameters.filestem = filestem;
filehandlingparameters.savelogfile = 1;
filehandlingparameters.showtext = 0;

EDres = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Present the results
    
ntotlength = size(EDres.irdirect,1);
allirtots = zeros(ntotlength,4);

allirtots(:,1) = EDres.irdirect + EDres.irgeom;
allirtots(:,2) = EDres.irdiff;
allirtots(:,3) = EDres.irhod{2};
allirtots(:,4) = EDres.irhod{3};

cumsumirtots = cumsum(allirtots.').';

tvec = 1/controlparameters.fs*[0:ntotlength-1].';

figure(1)
clf(1)
h = plot(tvec*1e3,allirtots(:,1)+allirtots(:,2),'-o');
set(h(1),'LineWidth',1,'MarkerSize',3);
g = get(h(1),'Parent');
set(g,'FontSize',14)
g = xlabel('Time   [ms]');
set(g,'FontSize',14)
g = ylabel('Impulse response   [-]');
set(g,'FontSize',14)
grid
axis([5 8 -0.1 0.1])
g = legend('Direct sound (w spec. refl.) and 1st-order diffr.');
set(g,'FontSize',14)

figure(2)
clf(2)
h = plot(tvec*1e3,allirtots(:,3),'-o');
set(h(1),'LineWidth',1,'MarkerSize',3);
g = get(h(1),'Parent');
set(g,'FontSize',14)
g = xlabel('Time   [ms]');
set(g,'FontSize',14)
g = ylabel('Impulse response   [-]');
set(g,'FontSize',14)
grid
axis([5 12 -0.001 0.003])
g = legend('Second-order diffr.');
set(g,'FontSize',14)

figure(3)
clf(3)
h = plot(tvec*1e3,allirtots(:,4),'-o');
set(h(1),'LineWidth',1,'MarkerSize',3);
g = get(h(1),'Parent');
set(g,'FontSize',14)
g = xlabel('Time   [ms]');
set(g,'FontSize',14)
g = ylabel('Impulse response   [-]');
set(g,'FontSize',14)
grid
axis([5 12 -0.0025 0.0005])
g = legend('Third-order diffr.');
set(g,'FontSize',14)


nfft = 4096;
fvec = controlparameters.fs/nfft*[0:nfft/2-1];
F = fft(cumsumirtots,nfft);

figure(4)
clf(4)
h = semilogx(fvec,20*log10(abs(F(1:nfft/2,:))));
for ii = 1:4
   set(h(ii),'LineWidth',2) 
end
g = get(h(1),'Parent');
set(g,'FontSize',14)
grid
g = xlabel('Frequency   [Hz]');
set(g,'FontSize',14)
g = ylabel('Frequency response magnitude   [dB]');
set(g,'FontSize',14)
g = legend('GA only','Incl. diffr.1','Incl. diffr.2','Incl. diffr.3');
axis([20 20000 -8 4])
