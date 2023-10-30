function [tftot_ESIEBEM,elapsedtimeESIEBEMpropagate,existingfilename] = ...
    EDpropagateESIEBEM(tftot_surface,Rinputdata,fieldpoints,...
    Sinputdata,directsound_fieldpoints,controlparameters,envdata,...
    EDversionnumber,filehandlingparameters)
% EDpropagateESIEBEM calculates the pressure at fieldpoints from the
% surface pressure.
%
% From v0.4 of the EDtoolbox, the input parameter list changed.
%
% Input parameters:
%   tftot_surface           [nsurfacepoints,nfreqs] The sound pressure at the surface points
%   Rinputdata              struct with fields
%                           .coordinates [nsurfacepoints,3]
%                           .nvecs       [nsurfacepoints,3]
%                           .weights     [nsurfacepoints,1]
%   fieldpoints             [nfieldpoints,3] The coordinates of the
%                           fieldpoints
%   Sinputdata
%   directsound_fieldpoints 0 or 1
%   controlparameters       struct, for which these fields are used:
%                           .frequencies
%                           .Rstart
%   envdata                 Input struct; field .cair is used here
%   EDversionnumber 
%   filehandlingparameters
% 
% Output parameters:
%   tftot_ESIEBEM           [nfrequencies,nreceivers,nsources]
%                           (if doaddsources = 0) or [nfrequencies,nreceivers]
%                           (if doaddsources = 1)
%   elapsedtimeESIEBEMpropagate
%                           This tells how long time was used inside this
%                           function. If an existing file was reused, then
%                           elapsedtimeESIEBEMpropagate has a second value which tells
%                           how much time was used for the existing file.
%   existingfilename        If an existing file was
%                           found that was reused, then the reused file
%                           name is given here. If no existing file could be 
%                           reused then this variable is empty. 
%
% Uses functions EDcalcdist, EDcalccosfi, EDrecycleresultfiles from EDtoolbox
% Uses function DataHash form Matlab Central
% 
% Peter Svensson 30 Oct. 2023 (peter.svensson@ntnu.no)
%
% [tftot_ESIEBEM,elapsedtimeESIEBEMpropagate,existingfilename] = ...
%    EDpropagateESIEBEM(tftot_surface,surfacepoints,surfacenvecs,fieldpoints,...
%    originalsourcecoordinates,directsound_fieldpoints,frequencies,envdata,...
%    EDversionnumber,filehandlingparameters)

% 29 Sep. 2023 Moved code from EDmain_convexESIEBEM into this function
% 3 Oct. 2023 Added the EDsettingnshash as input parameter
% 27 Oct. 2023 Added the fieldpoints to the datahash.

t00 = clock;

EDinputdatastruct = struct('tftot_surface',tftot_surface,...
    'fieldpoints',fieldpoints,'EDversionnumber',EDversionnumber);
EDinputdatahash = DataHash(EDinputdatastruct);

%---------------------------------------------------------------
% Sort out the file business: can an existing file be used?
% Then copy the existing file to a new copy. Should the data be saved in a file? 

if filehandlingparameters.suppressresultrecycling == 1
	foundmatch = 0;
	existingfilename = '';
else
	[foundmatch,existingfilename] = ... 
		EDrecycleresultfiles(filehandlingparameters.outputdirectory,...
		'_tfESIEBEM',EDinputdatahash);
end

desiredname = [filehandlingparameters.outputdirectory,filesep,...
	filehandlingparameters.filestem,'_tfESIEBEM.mat'];

if foundmatch == 1
	eval(['load(''',existingfilename,''')'])
	if ~strcmp(existingfilename,desiredname)
		copyfile(existingfilename,desiredname);
	end
	elapsedtimeESIEBEMpropagate_new = etime(clock,t00);
	elapsedtimeESIEBEMpropagate = [elapsedtimeESIEBEMpropagate_new elapsedtimeESIEBEMpropagate];
	return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the results for the surface sound pressure, and propagate to the
% field points = the original receiver points.
% 
% If directsound_fieldpoints = 1, then we also have to add the
% **free-field** component (not direct sound!), from the original source
% to the fieldpoints. 

surfacepoints = Rinputdata.coordinates;
surfacenvecs = Rinputdata.nvecs;

nfieldpoints = size(fieldpoints,1);
nsurfacepoints = size(surfacepoints,1);
nsources = size(Sinputdata.coordinates,1);
nfrequencies = length(controlparameters.frequencies);

[distances,cosfi] = EDcalccosfi(fieldpoints,surfacepoints,surfacenvecs);

if directsound_fieldpoints == 1
   ffdist = EDcalcdist(fieldpoints,Sinputdata.coordinates);
%    if any(Sinputdata.sourceamplitudes~=1)
%       error(['ERROR: sourceamplitudes have not been implemented yet in EDmain_convexESIEBEM']) 
%    end
end

if Sinputdata.doaddsources == 1 || nsources == 1
    tftot_ESIEBEM = zeros(nfrequencies,nfieldpoints);
else
    tftot_ESIEBEM = zeros(nfrequencies,nfieldpoints,nsources);    
end
kvec = 2*pi*controlparameters.frequencies(:)./envdata.cair;
kvec = kvec(:,ones(1,nsurfacepoints));
onesvec = ones(1,nfrequencies);

if Sinputdata.doaddsources == 1 || nsources == 1
    % tftot has size [nfreq,nsurfacerec]
    % tffieldpoints has size [nfreq,nfieldpoints]    
    %   becomes tftot
    for ii = 1:nfieldpoints    
        tfonefp = tftot_surface.*exp(-1i*kvec.*distances(:,ii*onesvec).')./distances(:,ii*onesvec).'.*(1i*kvec + 1./distances(:,ii*onesvec).').*cosfi(:,ii*onesvec).'.*Rinputdata.weights(:,onesvec).'/4/pi;
        tftot_ESIEBEM(:,ii) = sum(tfonefp.').';
    end
%    tftot_ESIEBEM = tftot_ESIEBEM;
    
    if directsound_fieldpoints == 1
        kvec = 2*pi*controlparameters.frequencies(:)./envdata.cair;
        kvec = kvec(:,ones(1,nfieldpoints));
        sourceamp = Sinputdata.sourceamplitudes.';
        if size(sourceamp,2) == 1
            sourceamp = sourceamp(:,ones(1,nsurfacepoints));
        end
        for jj = 1:nsources
            tftot_ESIEBEM = tftot_ESIEBEM + exp(-1i*kvec.*(ffdist(:,onesvec*jj).'-controlparameters.Rstart))./ffdist(:,onesvec*jj).'.*sourceamp(:,jj);
        end
    end
else
    % tftot has size [nfreq,nsurfacerec,nsources]
    % tffieldpoints has size [nfreq,nfieldpoints,nsources]    
    %   becomes tftot
    for ii = 1:nfieldpoints    
        transfermatrix = -exp(-1i*kvec.*distances(:,ii*onesvec).')./distances(:,ii*onesvec).'.*(1i*kvec + 1./distances(:,ii*onesvec).').*cosfi(:,ii*onesvec).'.*Rinputdata.weights(:,onesvec).'/4/pi;
        for jj = 1:nsources
            if Sinputdata.doallSRcombinations == 1 || jj == ii
                tfonefp = squeeze(tftot_ESIEBEM(:,:,jj)).*transfermatrix;
                tftot_ESIEBEM(:,ii,jj) = sum(tfonefp.').';
            end
        end
    end

    checkvalue = (size(tftot_ESIEBEM,2)==1) || (size(tftot_ESIEBEM,3)==1);
    if directsound_fieldpoints == 1
        for kk = 1:nfrequencies
            kvec = 2*pi*controlparameters.frequencies(kk*ones(nfieldpoints,nsources))./envdata.cair;
            if checkvalue == 1 
                tftot_ESIEBEM(kk,:,:) = squeeze(tftot_ESIEBEM(kk,:,:)).' + exp(-1i*kvec.*(ffdist-controlparameters.Rstart))./ffdist;            
            else
                tftot_ESIEBEM(kk,:,:) = squeeze(tftot_ESIEBEM(kk,:,:)) + exp(-1i*kvec.*(ffdist-controlparameters.Rstart))./ffdist;                            
            end
        end
    end
    
end

elapsedtimeESIEBEMpropagate = etime(clock,t00);

desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tfESIEBEM.mat'];

eval(['save(''',desiredname,''',''tftot_ESIEBEM'',''EDinputdatahash'',''elapsedtimeESIEBEMpropagate'');'])


