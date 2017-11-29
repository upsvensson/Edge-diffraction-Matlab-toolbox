function EDmain_convexESIE(geofiledata,Sindata,Rindata,envdata,controlparameters,filehandlingparameters)
% EDmain_convexESIE - Calculates the specular and edge diffraction IRs and/or TFs
% for convex, rigid scattering objects, using the image source method for
% the specular reflection, the explicit integral formulation for first-order
% diffraction, and the edge source integral equation for the second- and
% higher-order diffraction.
%
% Input parameters are six structs with fields as specified below:
%   geofiledata         .geoinputfile        (obligatory)
%                       .firstcornertoskip   (default: 1e6)
%   Sindata             .coordinates         (obligatory)
%                       .doaddsources        (default: 0 = no)
%   Rindata             .coordinates         (obligatory)
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
%                       .showtext        (default: 1)
%                       .savecadgeofile      (default: 0)
%                       .saveSRdatafiles     (default: 0)
%                       .saveeddatafile      (default: 0)
%                       .logfilename         (default: '')
%                       .lineending          (set automatically)
% 
% Peter Svensson 29 Nov. 2017 (peter.svensson@ntnu.no)
%
% EDmain_convex(geofiledata,Sindata,Rindata,envdata,controlparameters,filehandlingparameters);

% 24 Nov. 2017 A first version = a trimmed version of ESIE2main,
% specialized for the scattering from a convex object.
% 28 Nov. 2017 Removed the JJ. Added the handling of timing data and
% writing to a logfile.
% 29 Nov. 2017 Added the semicolon after fclose(fid). Made some adjustments to
% handle input matrices instead of a cad file.

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

[geofiledata,Sindata,Rindata,envdata,controlparameters,filehandlingparameters] = EDcheckinputstructs(geofiledata,Sindata,Rindata,envdata,controlparameters,filehandlingparameters,1);

if filehandlingparameters.savesetupfile == 1
    varlist = 'geofiledata Sindata Rindata envdata controlparameters filehandlingparameters';
    eval(['save ',filehandlingparameters.outputdirectory,filesep,'results',filesep,filehandlingparameters.filestem,'_setup.mat ',varlist])
end

if filehandlingparameters.showtext >= 1
	disp(' ');disp('####################################################################')
              disp('#  EDmain_convexESIE, v. 29 Nov. 2017')
end
if ~isempty(filehandlingparameters.logfilename)
    fid = fopen(filehandlingparameters.logfilename,'w');
    if fid == -1
        disp('This file is not possible to open - check that it isn''t opened by any program!')
    	return
    end
    fwrite(fid,[' ',filehandlingparameters.lineending],'char');   
    fwrite(fid,['####################################################################',filehandlingparameters.lineending],'char');
    fwrite(fid,['#  EDmain_convexESIE, v. 29 Nov. 2017',filehandlingparameters.lineending],'char');
    fwrite(fid,[' ',filehandlingparameters.lineending],'char');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the input CAD-file, or input matrices, and create the planedata struct

if isfield(geofiledata,'geoinputfile')
    if filehandlingparameters.showtext >= 1
        disp(' ');disp('   Creating the planedata struct from the CAD file')
    end

    t00 = clock;
    [planedata,extraCATTdata] = EDreadcad(geofiledata.geoinputfile,'circ',0);
    if isempty(strfind(planedata.modeltype,'convex_ext')) && isempty(strfind(planedata.modeltype,'singleplate'))
        error('ERROR: EDmain_convexESIE can only be used for convex scatterers, including a single thin plate')
    end
    if filehandlingparameters.savecadgeofile == 1
        desiredname = [filehandlingparameters.outputdirectory,filesep,'results',filesep,filehandlingparameters.filestem,'_cadgeo.mat'];
        eval(['save ',desiredname,' planedata extraCATTdata'])
    end
    t01 = etime(clock,t00);
    if ~isempty(filehandlingparameters.logfilename)
        fwrite(fid,['   EDreadcad, time: ',num2str(t01),' s',filehandlingparameters.lineending],'char');
    end
else
    if filehandlingparameters.showtext >= 1
        disp(' ');disp('   Creating the planedata struct from the input geometry matrices')
    end

    t00 = clock;
    planedata = EDreadgeomatrices(geofiledata.corners,geofiledata.planecorners);
    if isempty(strfind(planedata.modeltype,'convex_ext')) && isempty(strfind(planedata.modeltype,'singleplate'))
        error('ERROR: EDmain_convexESIE can only be used for convex scatterers, including a single thin plate')
    end
    t01 = etime(clock,t00);
    if ~isempty(filehandlingparameters.logfilename)
        fwrite(fid,['   EDreadgeomatrices, time: ',num2str(t01),' s',filehandlingparameters.lineending],'char');
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the planedata struct and create an edgedata struct

if filehandlingparameters.showtext >= 1	
	disp(' ');disp('   Creating the edgedata struct ')
end

t00 = clock;
[edgedata,planedata] = EDedgeo(planedata,geofiledata.firstcornertoskip,[],0,filehandlingparameters.showtext);
if filehandlingparameters.saveeddatafile == 1
    desiredname = [filehandlingparameters.outputdirectory,filesep,'results',filesep,filehandlingparameters.filestem,'_eddata.mat'];
    eval(['save ',desiredname,' planedata edgedata'])
end
t01 = etime(clock,t00);
if ~isempty(filehandlingparameters.logfilename)
    fwrite(fid,['   EDedgeo, time: ',num2str(t01),' s',filehandlingparameters.lineending],'char');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the Sdata and Rdata structs

if filehandlingparameters.showtext >= 1	
	disp(' ');disp('   Creating the Sdata struct ')
end

t00 = clock;
Sdata = EDSorRgeo(planedata,edgedata,Sindata.coordinates,'S',controlparameters.nedgepoints_visibility,filehandlingparameters.showtext);
if filehandlingparameters.saveSRdatafiles == 1
    desiredname = [filehandlingparameters.outputdirectory,filesep,'results',filesep,filehandlingparameters.filestem,'_Sdata.mat'];
    eval(['save ',desiredname,' Sdata'])    
end
t01 = etime(clock,t00);
if ~isempty(filehandlingparameters.logfilename)
    fwrite(fid,['   EDSorRgeo(S), time: ',num2str(t01),' s',filehandlingparameters.lineending],'char');
end

if filehandlingparameters.showtext >= 1	
	disp(' ');disp('   Creating the Rdata struct ')
end
t00 = clock;
Rdata = EDSorRgeo(planedata,edgedata,Rindata.coordinates,'R',controlparameters.nedgepoints_visibility,filehandlingparameters.showtext);
if filehandlingparameters.saveSRdatafiles == 1
    desiredname = [filehandlingparameters.outputdirectory,filesep,'results',filesep,filehandlingparameters.filestem,'_Rdata.mat'];
    eval(['save ',desiredname,' Rdata'])    
end
t01 = etime(clock,t00);
if ~isempty(filehandlingparameters.logfilename)
    fwrite(fid,['   EDSorRgeo(R), time: ',num2str(t01),' s',filehandlingparameters.lineending],'char');
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
if ~isempty(filehandlingparameters.logfilename)
    fwrite(fid,['   EDed2geo, time: ',num2str(t01),' s',filehandlingparameters.lineending],'char');
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
if ~isempty(filehandlingparameters.logfilename)
    fwrite(fid,['   EDinteg_submatrixstructure, time: ',num2str(t01),' s',filehandlingparameters.lineending],'char');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the HOD contribution with the integral equation

if filehandlingparameters.showtext >= 1	
	disp(' ');disp('   Calculating the HOD contribution with the FD integral equation')
end

t00 = clock;
[tfinteqdiff,timingdata,extraoutputdata] = EDintegralequation_convex_tf(filehandlingparameters,...
    envdata,planedata,edgedata,edgetoedgedata,Hsubmatrixdata,Sdata,Sindata.doaddsources,...
        Rdata,controlparameters);
desiredname = [filehandlingparameters.outputdirectory,filesep,'results',filesep,filehandlingparameters.filestem,'_tfinteq.mat'];
eval(['save ',desiredname,' tfinteqdiff extraoutputdata'])
t01 = etime(clock,t00);
if ~isempty(filehandlingparameters.logfilename)
    fwrite(fid,['   EDintegralequation_convex_tf, time: ',num2str(t01),' s. Parts, for one freq: ',num2str(timingdata),' s',filehandlingparameters.lineending],'char');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the paths for direct sound, first-order specular, and first-order
% diffraction paths.

if filehandlingparameters.showtext >= 1	
	disp(' ');disp('   Find the (first-order) GA paths, and tfs, one source at a time')
end

nsources = size(Sdata.sources,1);
nreceivers = size(Rdata.receivers,1);

for isou = 1:nsources
    if isou ==1
        t00 = clock;
    end
    [lengthNspecmatrix,lengthNdiffmatrix,singlediffcol,startindicessinglediff,endindicessinglediff,ndecimaldivider,PointertoIRcombs,IRoriginsfrom] = ...
    EDfindISEStree(planedata,edgedata,edgetoedgedata,Sdata.sources(isou,:),1,1,Sdata.visplanesfroms(:,isou),Sdata.vispartedgesfroms(:,isou),filehandlingparameters.showtext);
    if isou == 1
        t01 = etime(clock,t00);
        if ~isempty(filehandlingparameters.logfilename)
            fwrite(fid,['   EDfindISEStree, one source time: ',num2str(t01),' s',filehandlingparameters.lineending],'char');
        end
    end

    if isou ==1
        t00 = clock;
    end
    for irec = 1:nreceivers
        pathstruct = EDfindpathsISESx(planedata,edgedata,lengthNspecmatrix,lengthNdiffmatrix,singlediffcol,startindicessinglediff,...
            endindicessinglediff,ndecimaldivider,PointertoIRcombs,IRoriginsfrom,...
            Sdata.sources(isou,:),Rdata.receivers(irec,:),controlparameters.directsound,1,1,...
            controlparameters.nedgepoints_visibility,Sdata.visplanesfroms(:,isou),Rdata.visplanesfromr(:,irec),...
            Sdata.vispartedgesfroms(:,isou),Rdata.vispartedgesfromr(:,irec),filehandlingparameters.showtext);
            pathdata{isou,irec} = pathstruct;    
            
            [tfdirect_one,tfgeom_one,tfdiff_one] = EDmaketfs(envdata,planedata,edgedata,edgetoedgedata,...
                pathstruct,1,controlparameters.frequencies,controlparameters.Rstart,...                        
                [],controlparameters.directsound,filehandlingparameters.showtext);

            if Sindata.doaddsources == 1
                tfdirect(:,irec) = tfdirect(:,irec) + tfdirect_one; 
                tfgeom(:,irec)   = tfgeom(:,irec)   + tfgeom_one; 
                tfdiff(:,irec)   = tfdiff(:,irec)   + tfdiff_one; 
            elseif Sindata.doaddsources == 0
                tfdirect(:,irec,isou) =  tfdirect_one; 
                tfgeom(:,irec,isou)   = tfgeom_one; 
                tfdiff(:,irec,isou)   = tfdiff_one; 
            end
        
    end
    if isou == 1
        t01 = etime(clock,t00);
        if ~isempty(filehandlingparameters.logfilename)
            fwrite(fid,['   EDfindpathsISESx and EDmaketfs, one source, time: ',num2str(t01),' s',filehandlingparameters.lineending],'char');
        end
    end

end
desiredname = [filehandlingparameters.outputdirectory,filesep,'results',filesep,filehandlingparameters.filestem,'_tf.mat'];
eval(['save ',desiredname,' tfdirect tfgeom tfdiff'])

if ~isempty(filehandlingparameters.logfilename)
    fclose(fid);
end



