function EDmain_convexESIE_cornerspecial(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters)
% EDmain_convexESIE - Calculates the specular and edge diffraction IRs and/or TFs
% for convex, rigid scattering objects, using the image source method for
% the specular reflection, the explicit integral formulation for first-order
% diffraction, and the edge source integral equation for the second- and
% higher-order diffraction.
%
% Input parameters are six structs with fields as specified below:
%   geoinputdata         .geoinputfile        (obligatory)
%                       .firstcornertoskip   (default: 1e6)
%   Sinputdata             .coordinates         (obligatory)
%                       .doaddsources        (default: 0 = no)
%   Rinputdata             .coordinates         (obligatory)
%   envdata             .cair                (default: 344)
%                       .rhoair              (default: 1.21)
%   controlparameters   .fs                  (default: 44100)
%                       .directsound         (default: 1 = yes)
%                       .difforder           (default: 15)
%                       .docalctf            (default: 1)
%                       .docalcir            (default: 0)
%                       .Rstart              (default: 0)
%                       .frequencies         (obligatory, if docalcftf = 1)
%                       .discretizationtype  (default: 2 = G-L)
%                       .ngauss              (default: 16)
%   filehandlingparameters (optional)    
%                       .outputdirectory  (default: /result, in the folder of the geoinputfile)  
%                       .filestem        (default: name of the cad-file, with an underscore + running integer)
%                       .savesetupfile       (default: 1)
%                       .showtext            (default: 1)
%                       .savecadgeofile      (default: 0)
%                       .saveSRdatafiles     (default: 1)
%                       .saveeddatafile      (default: 1)
%                       .savepathsfile       (default: 0)
%                       .saveISEStree        (default: 0)
%                       .savelogfile         (default: 1)
%                       .savediff2result      (default: 0)
% 
% Peter Svensson 30 Nov. 2017 (peter.svensson@ntnu.no)
%
% EDmain_convex_cornerspecial(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters);

% 24 Nov. 2017 A first version = a trimmed version of ESIE2main,
% specialized for the scattering from a convex object.
% 28 Nov. 2017 Removed the JJ. Added the handling of timing data and
% writing to a logfile.
% 29 Nov. 2017 Added the semicolon after fclose(fid). Made some adjustments to
% handle input matrices instead of a cad file.
% 30 Nov. 2017 Fixed a missing part with calctfs and doaddsources = 1.
% Changed one field to savelogfile instead of logfilename. Removed the
% lineending field.

global POTENTIALISES ISCOORDS IVNDIFFMATRIX
global IVNSPECMATRIX ORIGINSFROM ISESVISIBILITY REFLORDER

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input data, assign default values if needed

if nargin < 6
    filehandlingparameters = struct('showtext',0);
    if nargin < 5
        error('ERROR: The first five input parameters to EDmain_convexESIE must be specified')  
    end
end

compstr = computer;
compstr = lower(compstr(1:3));
if compstr == 'mac'  
	lineending = 13;
elseif compstr == 'sun' || compstr == 'sol'            
	lineending = 10;
elseif compstr == 'pcw'
	lineending = [13,10];
else
    error('ERROR: Not implemented for this computer type yet')	
end

if filehandlingparameters.savelogfile == 1
    logfilename = [filehandlingparameters.outputdirectory,filesep,'results',filesep,filehandlingparameters.filestem,'_log.txt'];
end

[geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters] = EDcheckinputstructs(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters,1);

if filehandlingparameters.savesetupfile == 1
    varlist = 'geoinputdata Sinputdata Rinputdata envdata controlparameters filehandlingparameters';
    eval(['save ',filehandlingparameters.outputdirectory,filesep,'results',filesep,filehandlingparameters.filestem,'_setup.mat ',varlist])
end

if filehandlingparameters.showtext >= 1
	disp(' ');disp('####################################################################')
              disp('#  EDmain_convexESIE, v. 30 Nov. 2017')
end
if filehandlingparameters.savelogfile == 1
    fid = fopen(logfilename,'w');
    if fid == -1
        disp('This file is not possible to open - check that it isn''t opened by any program!')
    	return
    end
    fwrite(fid,[' ',lineending],'char');   
    fwrite(fid,['####################################################################',lineending],'char');
    fwrite(fid,['#  EDmain_convexESIE, v. 30 Nov. 2017',lineending],'char');
    fwrite(fid,[' ',lineending],'char');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the input CAD-file, or input matrices, and create the planedata struct

if isfield(geoinputdata,'geoinputfile')
    if filehandlingparameters.showtext >= 1
        disp(' ');disp('   Creating the planedata struct from the CAD file')
    end

    t00 = clock;
    [planedata,extraCATTdata] = EDreadcad(geoinputdata.geoinputfile,'circ',0);
    if isempty(strfind(planedata.modeltype,'convex_ext')) && isempty(strfind(planedata.modeltype,'singleplate'))
        error('ERROR: EDmain_convexESIE can only be used for convex scatterers, including a single thin plate')
    end
    if filehandlingparameters.savecadgeofile == 1
        desiredname = [filehandlingparameters.outputdirectory,filesep,'results',filesep,filehandlingparameters.filestem,'_cadgeo.mat'];
        eval(['save ',desiredname,' planedata extraCATTdata'])
    end
    t01 = etime(clock,t00);
    if filehandlingparameters.savelogfile == 1
        fwrite(fid,['   EDreadcad, time: ',num2str(t01),' s',lineending],'char');
    end
else
    if filehandlingparameters.showtext >= 1
        disp(' ');disp('   Creating the planedata struct from the input geometry matrices')
    end

    t00 = clock;
    planedata = EDreadgeomatrices(geoinputdata.corners,geoinputdata.planecorners);    
    t01 = etime(clock,t00);
    if filehandlingparameters.savelogfile == 1
        fwrite(fid,['   EDreadgeomatrices, time: ',num2str(t01),' s',lineending],'char');
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the planedata struct and create an edgedata struct

if filehandlingparameters.showtext >= 1	
	disp(' ');disp('   Creating the edgedata struct ')
end

t00 = clock;
[edgedata,planedata] = EDedgeo(planedata,geoinputdata.firstcornertoskip,[],0,filehandlingparameters.showtext);
if filehandlingparameters.saveeddatafile == 1
    desiredname = [filehandlingparameters.outputdirectory,filesep,'results',filesep,filehandlingparameters.filestem,'_eddata.mat'];
    eval(['save ',desiredname,' planedata edgedata'])
end
t01 = etime(clock,t00);
if filehandlingparameters.savelogfile == 1
    fwrite(fid,['   EDedgeo, time: ',num2str(t01),' s',lineending],'char');
end
if isempty(strfind(planedata.modeltype,'convex_ext')) && isempty(strfind(planedata.modeltype,'singleplate'))
    error('ERROR: EDmain_convexESIE can only be used for convex scatterers, including a single thin plate')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the Sdata and Rdata structs

if filehandlingparameters.showtext >= 1	
	disp(' ');disp('   Creating the Sdata struct ')
end

t00 = clock;
Sdata = EDSorRgeo(planedata,edgedata,Sinputdata.coordinates,'S',controlparameters.nedgepoints_visibility,filehandlingparameters.showtext);
if filehandlingparameters.saveSRdatafiles == 1
    desiredname = [filehandlingparameters.outputdirectory,filesep,'results',filesep,filehandlingparameters.filestem,'_Sdata.mat'];
    eval(['save ',desiredname,' Sdata'])    
end
t01 = etime(clock,t00);
if filehandlingparameters.savelogfile == 1
    fwrite(fid,['   EDSorRgeo(S), time: ',num2str(t01),' s',lineending],'char');
end

if filehandlingparameters.showtext >= 1	
	disp(' ');disp('   Creating the Rdata struct ')
end
t00 = clock;
Rdata = EDSorRgeo(planedata,edgedata,Rinputdata.coordinates,'R',controlparameters.nedgepoints_visibility,filehandlingparameters.showtext);
if filehandlingparameters.saveSRdatafiles == 1
    desiredname = [filehandlingparameters.outputdirectory,filesep,'results',filesep,filehandlingparameters.filestem,'_Rdata.mat'];
    eval(['save ',desiredname,' Rdata'])    
end
t01 = etime(clock,t00);
if filehandlingparameters.savelogfile == 1
    fwrite(fid,['   EDSorRgeo(R), time: ',num2str(t01),' s',lineending],'char');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add edge-to-edge visibility data

if filehandlingparameters.showtext >= 1	
	disp(' ');disp('   Creating the edgetoedgedata struct ')
end
t00 = clock;
edgetoedgedata = EDed2geo(edgedata,planedata,Sdata,Rdata,1,2,controlparameters.nedgepoints_visibility,filehandlingparameters.showtext);    
if filehandlingparameters.saveeddatafile == 1
    desiredname = [filehandlingparameters.outputdirectory,filesep,'results',filesep,filehandlingparameters.filestem,'_eddata.mat'];
    eval(['save ',desiredname,' planedata edgedata edgetoedgedata'])
end
t01 = etime(clock,t00);
if filehandlingparameters.savelogfile == 1
    fwrite(fid,['   EDed2geo, time: ',num2str(t01),' s',lineending],'char');
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare for the integral equation: set up the submatrix structure

if filehandlingparameters.showtext >= 1	
	disp(' ');disp('   Creating the integral equation Hsubmatrixdata struct ')
end

t00 = clock;
Hsubmatrixdata = EDinteg_submatrixstructure(edgedata.edgelengthvec,edgedata.closwedangvec,...
    controlparameters.ngauss,controlparameters.discretizationtype,edgetoedgedata,edgedata.planesatedge,filehandlingparameters.showtext);
if filehandlingparameters.savesubmatrixdata == 1
    desiredname = [filehandlingparameters.outputdirectory,filesep,'results',filesep,filehandlingparameters.filestem,'_submatrixdata.mat'];
    eval(['save ',desiredname,' Hsubmatrixdata'])
end
t01 = etime(clock,t00);
if filehandlingparameters.savelogfile == 1
    fwrite(fid,['   EDinteg_submatrixstructure, time: ',num2str(t01),' s',lineending],'char');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the HOD contribution with the integral equation

if filehandlingparameters.showtext >= 1	
	disp(' ');disp('   Calculating the HOD contribution with the FD integral equation')
end

t00 = clock;
[tfinteqdiff,timingdata,extraoutputdata] = EDintegralequation_convex_tf(filehandlingparameters,...
    envdata,planedata,edgedata,edgetoedgedata,Hsubmatrixdata,Sdata,Sinputdata.doaddsources,...
        Rdata,controlparameters);
desiredname = [filehandlingparameters.outputdirectory,filesep,'results',filesep,filehandlingparameters.filestem,'_tfinteq.mat'];
eval(['save ',desiredname,' tfinteqdiff extraoutputdata'])
t01 = etime(clock,t00);
if filehandlingparameters.savelogfile == 1
    fwrite(fid,['   EDintegralequation_convex_tf, time: ',num2str(t01),' s. Parts, for one freq: ',num2str(timingdata),' s',lineending],'char');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the paths for direct sound, first-order specular, and first-order
% diffraction paths.

if filehandlingparameters.showtext >= 1	
	disp(' ');disp('   Find the (first-order) GA paths, and tfs, one source at a time')
end

nsources = size(Sdata.sources,1);
nreceivers = size(Rdata.receivers,1);
nfrequencies = length(controlparameters.frequencies);

if controlparameters.docalctf == 1
    if Sinputdata.doaddsources == 1
       tfdirect = zeros(nfrequencies,nreceivers);
       tfgeom = zeros(nfrequencies,nreceivers);
       tfdiff = zeros(nfrequencies,nreceivers);        
    elseif Sinputdata.doaddsources == 0
       tfdirect = zeros(nfrequencies,nreceivers,nsources);
       tfgeom = zeros(nfrequencies,nreceivers,nsources);
       tfdiff = zeros(nfrequencies,nreceivers,nsources);        
    end
end
if controlparameters.docalcir == 1
    if Sinputdata.doaddsources == 1
       irdirect = zeros(ntimesteps,nreceivers);
       irgeom = zeros(ntimesteps,nreceivers);
       irdiff = zeros(ntimesteps,nreceivers);
    elseif Sinputdata.doaddsources == 0
       irdirect = zeros(ntimesteps,nreceivers,nsources);
       irgeom = zeros(ntimesteps,nreceivers,nsources);
       irdiff = zeros(ntimesteps,nreceivers,nsources);
    end
end

for isou = 1:nsources
    if isou ==1
        t00 = clock;
    end
%     [lengthNspecmatrix,lengthNdiffmatrix,singlediffcol,startindicessinglediff,endindicessinglediff,ndecimaldivider,PointertoIRcombs,IRoriginsfrom] = ...
%     EDfindISEStree(planedata,edgedata,edgetoedgedata,Sdata.sources(isou,:),1,1,Sdata.visplanesfroms(:,isou),Sdata.vispartedgesfroms(:,isou),filehandlingparameters.showtext);
    [lengthNspecmatrix,lengthNdiffmatrix,singlediffcol,startindicessinglediff,endindicessinglediff,ndecimaldivider,PointertoIRcombs,IRoriginsfrom] = ...
    EDfindISEStree(planedata,edgedata,edgetoedgedata,Sdata.sources(isou,:),2,2,Sdata.visplanesfroms(:,isou),Sdata.vispartedgesfroms(:,isou),filehandlingparameters.showtext);
disp('WARNING: temp. changed difforder to 2; should be 1')
    if filehandlingparameters.saveISEStree == 1
        desiredname = [filehandlingparameters.outputdirectory,filesep,'results',filesep,filehandlingparameters.filestem,'_',int2str(isou),'_ISEStree.mat'];
        eval(['save ',desiredname,' POTENTIALISES ISCOORDS IVNDIFFMATRIX IVNSPECMATRIX ORIGINSFROM ISESVISIBILITY REFLORDER ', ...
        ' lengthNspecmatrix lengthNdiffmatrix singlediffcol startindicessinglediff endindicessinglediff ndecimaldivider PointertoIRcombs IRoriginsfrom'])                 
    end
    if isou == 1
        t01 = etime(clock,t00);
        if filehandlingparameters.savelogfile == 1
            fwrite(fid,['   EDfindISEStree, one source, time: ',num2str(t01),' s',lineending],'char');
        end
    end

    if isou ==1
        t00 = clock;
    end
    for irec = 1:nreceivers
%         pathstruct = EDfindpathsISESx(planedata,edgedata,lengthNspecmatrix,lengthNdiffmatrix,singlediffcol,startindicessinglediff,...
%             endindicessinglediff,ndecimaldivider,PointertoIRcombs,IRoriginsfrom,...
%             Sdata.sources(isou,:),Rdata.receivers(irec,:),controlparameters.directsound,1,1,...
%             controlparameters.nedgepoints_visibility,Sdata.visplanesfroms(:,isou),Rdata.visplanesfromr(:,irec),...
%             Sdata.vispartedgesfroms(:,isou),Rdata.vispartedgesfromr(:,irec),filehandlingparameters.showtext);
%             pathdata{isou,irec} = pathstruct;    
disp('WARNING: temp. changed difforder to 2; should be 1')
        pathstruct = EDfindpathsISESx(planedata,edgedata,lengthNspecmatrix,lengthNdiffmatrix,singlediffcol,startindicessinglediff,...
            endindicessinglediff,ndecimaldivider,PointertoIRcombs,IRoriginsfrom,...
            Sdata.sources(isou,:),Rdata.receivers(irec,:),controlparameters.directsound,2,2,...
            controlparameters.nedgepoints_visibility,Sdata.visplanesfroms(:,isou),Rdata.visplanesfromr(:,irec),...
            Sdata.vispartedgesfroms(:,isou),Rdata.vispartedgesfromr(:,irec),filehandlingparameters.showtext);
            pathdata{isou,irec} = pathstruct;   
            
            if filehandlingparameters.savepathsfile == 1
                desiredname = [filehandlingparameters.outputdirectory,filesep,'results',filesep,filehandlingparameters.filestem,'_',int2str(isou),'_',int2str(irec),'_edpaths.mat'];
                eval(['save ',desiredname,' pathstruct'])                 
            end
            
            if controlparameters.docalctf == 1
%                 [tfdirect_one,tfgeom_one,tfdiff_one] = EDmaketfs(envdata,planedata,edgedata,edgetoedgedata,...
%                     pathstruct,1,controlparameters.frequencies,controlparameters.Rstart,...                        
%                     [],controlparameters.directsound,filehandlingparameters.showtext);
disp('WARNING: temp. changed specorder to 2; should be 1')
                [tfdirect_one,tfgeom_one,tfdiff_one] = EDmaketfs_cornerspecial(envdata,planedata,edgedata,edgetoedgedata,...
                    pathstruct,2,controlparameters.frequencies,controlparameters.Rstart,...                        
                    [],controlparameters.directsound,filehandlingparameters.showtext);

                if Sinputdata.doaddsources == 1
                    tfdirect(:,irec) = tfdirect(:,irec) + tfdirect_one; 
                    tfgeom(:,irec)   = tfgeom(:,irec)   + tfgeom_one; 
                    tfdiff(:,irec)   = tfdiff(:,irec)   + tfdiff_one; 
                elseif Sinputdata.doaddsources == 0
                    tfdirect(:,irec,isou) =  tfdirect_one; 
                    tfgeom(:,irec,isou)   = tfgeom_one; 
                    tfdiff(:,irec,isou)   = tfdiff_one; 
                end
            end
            if controlparameters.docalcir == 1
                disp('WARNING: calcirs not implemented yet');
            end
            
    end
    if isou == 1
        t01 = etime(clock,t00);
        if filehandlingparameters.savelogfile == 1
            fwrite(fid,['   EDfindpathsISESx and EDmaketfs, one source, time: ',num2str(t01),' s',lineending],'char');
        end
    end

end
desiredname = [filehandlingparameters.outputdirectory,filesep,'results',filesep,filehandlingparameters.filestem,'_tf.mat'];
eval(['save ',desiredname,' tfdirect tfgeom tfdiff'])

if filehandlingparameters.savelogfile == 1
    fclose(fid);
end



