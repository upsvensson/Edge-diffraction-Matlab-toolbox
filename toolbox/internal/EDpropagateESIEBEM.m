function [tftot_ESIEBEM,elapsedtimeESIEBEMpropagate,existingfilename] = ...
    EDpropagateESIEBEM(tftot_surface,Rinputdata,fieldpoints,...
    Sinputdata,directsound_fieldpoints,controlparameters,envdata,...
    EDversionnumber,filehandlingparameters)
% EDpropagateESIEBEM - Computes the sound pressure at field points from the surface pressure (ESIEBEM).
% The sound pressure at the surface points (computed by the ESIE) is
% propagated to the field points with the Kirchhoff-Helmholtz integral.
% If directsound_fieldpoints = 1, the free-field sound from the source is
% added. Only a single source, or added sources, can be handled. The
% results are saved in the stage file <filestem>_tfESIEBEM.mat, which is
% reused when the input data are the same.
%
% Input parameters:
%   tftot_surface   [nfrequencies,nsurfacepoints] the sound pressure at the
%                   surface points
%   Rinputdata      struct with the surface points: .coordinates
%                   [nsurfacepoints,3], .nvecs [nsurfacepoints,3] and
%                   .weights [nsurfacepoints,1] (from EDgensurfreceivers)
%   fieldpoints     [nfieldpoints,3] the coordinates of the field points
%                   (the receivers given by the user)
%   Sinputdata      struct; the fields .coordinates, .sourceamplitudes,
%                   .doaddsources and .doallSRcombinations are used.
%   directsound_fieldpoints  0 or 1
%   controlparameters  struct; the fields .frequencies and .Rstart are
%                   used.
%   envdata         struct; the field .cair is used.
%   EDversionnumber  The version number of the EDtoolbox (stored in the
%                   stage file).
%   filehandlingparameters
%                   struct; the fields .outputdirectory, .filestem and
%                   .suppressresultrecycling are used.
%
% Output parameters:
%   tftot_ESIEBEM   [nfrequencies,nfieldpoints] the sound pressure at the
%                   field points
%   elapsedtimeESIEBEMpropagate
%                   The time used by this function, in seconds. If an
%                   existing stage file was reused, a second value gives
%                   the time it took to compute that file.
%   existingfilename  The name of the stage file that was reused, or empty
%                   if the data were computed.
%
% Uses the functions EDcalccosfi, EDcalcdist, EDrecycleresultfiles,
% EDstageversions from the EDtoolbox
% Uses the function DataHash from Matlab Central
%
% Peter Svensson (peter.svensson@ntnu.no) 14 Oct. 2026
%
% [tftot_ESIEBEM,elapsedtimeESIEBEMpropagate,existingfilename] = ...
%    EDpropagateESIEBEM(tftot_surface,Rinputdata,fieldpoints,Sinputdata,...
%    directsound_fieldpoints,controlparameters,envdata,EDversionnumber,...
%    filehandlingparameters);

% 29 Sep. 2023 Moved code from EDmain_convexESIEBEM into this function
% 3 Oct. 2023 Added the EDsettingnshash as input parameter
% 27 Oct. 2023 Added the fieldpoints to the datahash.
% 14 Oct. 2026 The EDinputdatahash was incomplete: it contained only the
% surface pressure and the field points. Now it also contains the surface
% points, normal vectors and weights, the source data, the direct sound
% setting, the frequencies, Rstart and cair, and a stageversion from
% EDstageversions (instead of the overall EDversionnumber). Before, a
% change of only controlparameters.directsound led to the reuse of an old
% _tfESIEBEM file. The branch for several, not added, sources was wrong
% (it used the output matrix instead of tftot_surface, had a different
% sign, and ignored the source amplitudes). It is not reachable from
% EDmain_convex, and has been replaced by an error message.
% load/save without eval.
% 14 Oct. 2026 The header was updated to the standard structure (by Claude,
% prompted by PS).

t00 = clock;

% The version of this function, and of the subroutines it uses, are
% defined in EDstageversions. They are part of the hash, so that saved
% stage files are not reused after a change in the computations.
stageversion = EDstageversions('EDpropagateESIEBEM');

EDinputdatastruct = struct('stage','EDpropagateESIEBEM','stageversion',stageversion,...
    'tftot_surface',tftot_surface,'fieldpoints',fieldpoints,...
    'surfacepoints',Rinputdata.coordinates,'surfacenvecs',Rinputdata.nvecs,...
    'surfaceweights',Rinputdata.weights,...
    'sources',Sinputdata.coordinates,'sourceamplitudes',Sinputdata.sourceamplitudes,...
    'doaddsources',Sinputdata.doaddsources,'doallSRcombinations',Sinputdata.doallSRcombinations,...
    'directsound_fieldpoints',directsound_fieldpoints,...
    'frequencies',controlparameters.frequencies,'Rstart',controlparameters.Rstart,...
    'cair',envdata.cair);
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
	filedata = load(existingfilename,'tftot_ESIEBEM','elapsedtimeESIEBEMpropagate');
	tftot_ESIEBEM = filedata.tftot_ESIEBEM;
	elapsedtimeESIEBEMpropagate = filedata.elapsedtimeESIEBEMpropagate;
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

if ~(Sinputdata.doaddsources == 1 || nsources == 1)
    error('ERROR: EDpropagateESIEBEM has not been implemented for several sources that are not added (doaddsources = 0)')
end

[distances,cosfi] = EDcalccosfi(fieldpoints,surfacepoints,surfacenvecs);

if directsound_fieldpoints == 1
   ffdist = EDcalcdist(fieldpoints,Sinputdata.coordinates);
%    if any(Sinputdata.sourceamplitudes~=1)
%       error(['ERROR: sourceamplitudes have not been implemented yet in EDmain_convexESIEBEM']) 
%    end
end

tftot_ESIEBEM = zeros(nfrequencies,nfieldpoints);
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
end

elapsedtimeESIEBEMpropagate = etime(clock,t00);

desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tfESIEBEM.mat'];

save(desiredname,'tftot_ESIEBEM','EDinputdatahash','elapsedtimeESIEBEMpropagate','stageversion','EDversionnumber');


