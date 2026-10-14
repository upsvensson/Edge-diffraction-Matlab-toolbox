function [irhod,elapsedtimemakeirhod,existingfilename] = EDmakeHODirs(hodpaths,...
    hodpathsalongplane,edgedata,edgetoedgedata,Sdata,Rdata,envdata,...
    controlparameters,EDversionnumber,filehandlingparameters)
% EDmakeHODirs - Computes the second- and higher-order diffraction IRs.
% Computes the impulse responses for the paths in hodpaths, with the
% discretized formulation of EDwedge2nd (second order) and EDwedgeN (third
% order and higher), and sums them, or keeps the diffraction orders
% separate. The results are saved in the stage file <filestem>_irhod.mat,
% which is reused when the input data are the same.
%
% Input parameters:
%   hodpaths, hodpathsalongplane  from EDfindHODpaths
%   edgedata, edgetoedgedata, Sdata, Rdata  structs, see help EDdatastructs
%   envdata         struct; the field .cair is used.
%   controlparameters  struct; the fields .difforder, .fs, .Rstart,
%                   .HODelemsize (the edge element size for each
%                   diffraction order) and .savealldifforders are used.
%   EDversionnumber  The version number of the EDtoolbox (stored in the
%                   stage file).
%   filehandlingparameters
%                   struct; the fields .outputdirectory, .filestem, .showtext and
%                   .suppressresultrecycling are used.
%
% Output parameters:
%   irhod           Cell array. If controlparameters.savealldifforders = 0,
%                   irhod{1} has the sum of all orders; if it is 1,
%                   irhod{N} has diffraction order N. Each cell has the size
%                   [nsamples,nreceivers,nsources], or [nsamples,nreceivers]
%                   if the sources are added (Sdata.doaddsources = 1).
%   elapsedtimemakeirhod
%                   The time used by this function, in seconds. If an
%                   existing stage file was reused, a second value gives
%                   the time it took to compute that file.
%   existingfilename  The name of the stage file that was reused, or empty
%                   if the data were computed.
%
% See also EDdatastructs, EDstageversions
%
% Peter Svensson (peter.svensson@ntnu.no) 15 Oct. 2026
%
% [irhod,elapsedtimemakeirhod,existingfilename] = EDmakeHODirs(hodpaths,...
%    hodpathsalongplane,edgedata,edgetoedgedata,Sdata,Rdata,envdata,...
%    controlparameters,EDversionnumber,filehandlingparameters);
    
% 8 Dec. 2006 First version
% 1 May 2017 Fixed a bug which gave an erroneous boosting of some
%            higher-order diffraction for thin planes.
% 1 May 2017 Fixed a bug which gave an error when there were specular
% reflections between diffractions. At least fixed for second-order
% diffraction. Possible error still: is the boosting done correctly for
% such cases?
% 12 Feb 2018 Adapted to the new output parameter of EDB1wedge1st_int.
% 16 Feb 2018 Condensed version from EDB1makeirs
% 16 Feb 2018 Renamed to irhod instead of hodir. Introduced the
% doaddsources and sourceamplitudes input parameters. Changed the calls
% from EDB1 functions to EDfunctions.
% 21 Mar 2018 Introduced new input parameter: savealldifforders
% 29 May 2018 Fixed small bug around line 300 for thin planes.
% 21 May 2019 Added the new input parameter hodpathsalongplane and modified
% the calls to EDwedge2nd and EDwedgeN accordingly.
% 24 Sep. 2022 Changed how the definition of the input data struct was
% defined. The previous version created a huge struct and failed to
% identify existing files that could be reused.
% 29 Sep. 2023 Implemented version 2 of this function while maintaining
% compatibility with the old "version 1". v2 moves the check if an existing
% file can be reused inside this function. Also updated load and save to
% the function call form, which avoids problems with spaces in file names.
% 27 Oct. 2023 Substantial change to the input parameter list.
% 29 Oct. 2023 Some more change to the input parameter list. elemsize is
% now a field in controlparameters.
% 29 Nov. 2023 Fixed small bug: the parameter cair was changed to
% envdata.cair for difforders 3 and higher.
% 12 Oct. 2026 New EDinputdatahash: chained to the hashes of edgedata,
% edgetoedgedata, Sdata and Rdata (instead of the complete structs), with a
% local stageversion instead of the overall EDversionnumber. Since Sdata is
% no longer hashed as a whole, sourceamplitudes and doaddsources are hashed
% explicitly; of envdata, only cair is used. load/save without eval.
% 14 Oct. 2026 The stageversion is taken from EDstageversions, and it
% also contains the versions of the subroutines that this function uses.
% 14 Oct. 2026 Fixed a bug: when the sources were not added (doaddsources
% = 0), the source amplitudes were not applied. Version 2 of this function
% in EDstageversions.
% 14 Oct. 2026 The header was updated to the standard structure (by Claude,
% prompted by PS).
% 15 Oct. 2026 The list of called functions was replaced by a See also line.

t00 = clock;

showtext = filehandlingparameters.showtext;
difforder = controlparameters.difforder;
fs = controlparameters.fs;
Rstart = controlparameters.Rstart;
elemsize = controlparameters.HODelemsize;

% The version of this function, and of the subroutines it uses, are
% defined in EDstageversions. They are part of the hash, so that saved
% stage files are not reused after a change in the computations.
stageversion = EDstageversions('EDmakeHODirs');

EDinputdatastruct = struct('stage','EDmakeHODirs','stageversion',stageversion,...
    'edgedatahash',EDparenthash(edgedata),'ed2datahash',EDparenthash(edgetoedgedata),...
    'Sdatahash',EDparenthash(Sdata),'Rdatahash',EDparenthash(Rdata),...
    'difforder',difforder,'fs',fs,'Rstart',Rstart,'cair',envdata.cair,...
    'savealldifforders',controlparameters.savealldifforders,...
    'sourceamplitudes',Sdata.sourceamplitudes,'doaddsources',Sdata.doaddsources);
% hodpaths and hodpathsalongplane are cell variables, so they are added as
% fields (struct(...) with a cell value would create a struct array)
EDinputdatastruct.hodpaths = hodpaths;
EDinputdatastruct.hodpathsalongplane = hodpathsalongplane;
EDinputdatastruct.elemsize = elemsize;
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
		'_irhod',EDinputdatahash);
end

desiredname = [filehandlingparameters.outputdirectory,filesep,...
	filehandlingparameters.filestem,'_irhod.mat'];

if foundmatch == 1
	filedata = load(existingfilename,'irhod','elapsedtimemakeirhod');
	irhod = filedata.irhod;
	elapsedtimemakeirhod = filedata.elapsedtimemakeirhod;
	if ~strcmp(existingfilename,desiredname)
		copyfile(existingfilename,desiredname);
	end
	elapsedtimemakeirhod_new = etime(clock,t00);
	elapsedtimemakeirhod = [elapsedtimemakeirhod_new elapsedtimemakeirhod];
	return
end

%----------------------------------------------------------------
% No existing file can be used

global BIGEDGESTEPMATRIX 

nyvec = pi./(2*pi - edgedata.closwedangvec);    

reftoshortlistE = edgetoedgedata.reftoshortlistE;

nsources = size(Sdata.coordinates,1);
nreceivers = size(Rdata.coordinates,1);

%--------------------------------------------------------------------------

firstcomponentdone = 0;
cellcounter = 1;

for isou = 1:nsources
    for irec = 1:nreceivers

        S = Sdata.coordinates(isou,:);
        R = Rdata.coordinates(irec,:);

        for Ndifforder = 2:difforder
            if showtext >= 2
                disp(['   Diffraction order ',int2str(Ndifforder)])    
            end

            if ~isempty(hodpaths{Ndifforder})
                
%                 disp(['   Paths found for order ',int2str(Ndifforder)])
                if controlparameters.savealldifforders == 1
                   cellcounter = Ndifforder; 
                end

                % 21 May 2019 Used the new input parameter
                % hodpathsalongplane
%                 pathalongplane = ones(1,Ndifforder-1);
                pathalongplane = hodpathsalongplane{Ndifforder,irec,isou};
                bc = ones(1,Ndifforder);

                % Calculate some general parameters that are shared for all
                % N-diffraction calculations

                divmin = envdata.cair/(fs*elemsize(Ndifforder));
                ndivvec = ceil(abs( edgedata.edgelengthvec.' )/divmin);
                dzvec = (edgedata.edgelengthvec.')./ndivvec;

                ncylrows = 4*(Ndifforder-1);

                %------------------------------------------------------------------
                % edgepattern will be a matrix which contains
                % the 2,3,4,... edge numbers for each path

                edgepatternlist = hodpaths{Ndifforder,irec,isou}; 
                ncomponents = size(edgepatternlist,1);                

                lastndivcomb = zeros(1,Ndifforder);

                for ii = 1:ncomponents

                    if showtext >= 4
                            disp(['      Combination no. ',int2str(ii),' of ',int2str(ncomponents),': ',int2str(edgepatternlist(ii,:))]) 
                    elseif showtext == 3
                        if round(ii/ceil(ncomponents/10))*ceil(ncomponents/10) == ii
                            disp(['      Combination no. ',int2str(ii),' of ',int2str(ncomponents)]) 
                        end
                    end
                    
                    if controlparameters.savealldifforders == 1
                        if ii == 1 && isou == 1 && irec == 1
                            firstcomponentdone = 0; 
                        end            
                    end

                    edgepattern = edgepatternlist(ii,:);

                    if Ndifforder >= 2

                        newndivvec = ndivvec(edgepattern);

                        if any(lastndivcomb~=newndivvec)
                            if Ndifforder == 2
                                ivmatrix = EDcreindexmatrix(newndivvec);
                            else
                                ivmatrix = EDcreindexmatrix(newndivvec(2:end));                        
                            end
                            nedgeelcombs = size(ivmatrix,1);
                            if Ndifforder == 2
                                BIGEDGESTEPMATRIX = (double(ivmatrix)-0.5)./newndivvec(uint8(ones(nedgeelcombs,1)),:);
                            else
                                BIGEDGESTEPMATRIX = (double(ivmatrix)-0.5)./newndivvec(uint8(ones(nedgeelcombs,1)),2:end);
                            end
                            
                            clear ivmatrix
                            lastndivcomb = newndivvec;

                        end
                    end  

                    firstedge = edgepattern(1);
                    lastedge = edgepattern(end);

                    cylS = [Sdata.rSsho(Sdata.reftoshortlistS(firstedge,isou)),...
                        Sdata.thetaSsho(Sdata.reftoshortlistS(firstedge,isou)),...
                        Sdata.zSsho(Sdata.reftoshortlistS(firstedge,isou))];

                    cylR = [Rdata.rRsho(Rdata.reftoshortlistR(lastedge,irec)),...
                        Rdata.thetaRsho(Rdata.reftoshortlistR(lastedge,irec)),...
                        Rdata.zRsho(Rdata.reftoshortlistR(lastedge,irec))];

                        % Pick out the edge-to-edge coordinates for this specific
                        % edge pair/N-let.                        

                    if Ndifforder >= 2
                        index1 = reftoshortlistE(edgepattern(2),edgepattern(1));
                        cylE2_r1 = [edgetoedgedata.re1sho(index1,:) edgetoedgedata.thetae1sho(index1,:) edgetoedgedata.ze1sho(index1,:);edgetoedgedata.re2sho(index1,:) edgetoedgedata.thetae2sho(index1,:) edgetoedgedata.ze2sho(index1,:)];
                        index2 = reftoshortlistE(edgepattern(1),edgepattern(2));
                        cylE1_r2 = [edgetoedgedata.re1sho(index2,:) edgetoedgedata.thetae1sho(index2,:) edgetoedgedata.ze1sho(index2,:);edgetoedgedata.re2sho(index2,:) edgetoedgedata.thetae2sho(index2,:) edgetoedgedata.ze2sho(index2,:)];
                    end
                    if Ndifforder >= 3
                        index1 = reftoshortlistE(edgepattern(3),edgepattern(2));
                        cylE3_r2 = [edgetoedgedata.re1sho(index1,:) edgetoedgedata.thetae1sho(index1,:) edgetoedgedata.ze1sho(index1,:);edgetoedgedata.re2sho(index1,:) edgetoedgedata.thetae2sho(index1,:) edgetoedgedata.ze2sho(index1,:)];
                        index2 = reftoshortlistE(edgepattern(2),edgepattern(3));
                        cylE2_r3 = [edgetoedgedata.re1sho(index2,:) edgetoedgedata.thetae1sho(index2,:) edgetoedgedata.ze1sho(index2,:);edgetoedgedata.re2sho(index2,:) edgetoedgedata.thetae2sho(index2,:) edgetoedgedata.ze2sho(index2,:)];
                    end
                    if Ndifforder >= 4
                        index1 = reftoshortlistE(edgepattern(4),edgepattern(3));
                        cylE4_r3 = [edgetoedgedata.re1sho(index1,:) edgetoedgedata.thetae1sho(index1,:) edgetoedgedata.ze1sho(index1,:);edgetoedgedata.re2sho(index1,:) edgetoedgedata.thetae2sho(index1,:) edgetoedgedata.ze2sho(index1,:)];
                        index2 = reftoshortlistE(edgepattern(3),edgepattern(4));
                        cylE3_r4 = [edgetoedgedata.re1sho(index2,:) edgetoedgedata.thetae1sho(index2,:) edgetoedgedata.ze1sho(index2,:);edgetoedgedata.re2sho(index2,:) edgetoedgedata.thetae2sho(index2,:) edgetoedgedata.ze2sho(index2,:)];
                    end
                    if Ndifforder >= 5
                        index1 = reftoshortlistE(edgepattern(5),edgepattern(4));
                        cylE5_r4 = [edgetoedgedata.re1sho(index1,:) edgetoedgedata.thetae1sho(index1,:) edgetoedgedata.ze1sho(index1,:);edgetoedgedata.re2sho(index1,:) edgetoedgedata.thetae2sho(index1,:) edgetoedgedata.ze2sho(index1,:)];
                        index2 = reftoshortlistE(edgepattern(4),edgepattern(5));
                        cylE4_r5 = [edgetoedgedata.re1sho(index2,:) edgetoedgedata.thetae1sho(index2,:) edgetoedgedata.ze1sho(index2,:);edgetoedgedata.re2sho(index2,:) edgetoedgedata.thetae2sho(index2,:) edgetoedgedata.ze2sho(index2,:)];
                    end
                    if Ndifforder >= 6
                        index1 = reftoshortlistE(edgepattern(6),edgepattern(5));
                        cylE6_r5 = [edgetoedgedata.re1sho(index1,:) edgetoedgedata.thetae1sho(index1,:) edgetoedgedata.ze1sho(index1,:);edgetoedgedata.re2sho(index1,:) edgetoedgedata.thetae2sho(index1,:) edgetoedgedata.ze2sho(index1,:)];
                        index2 = reftoshortlistE(edgepattern(5),edgepattern(6));
                        cylE5_r6 = [edgetoedgedata.re1sho(index2,:) edgetoedgedata.thetae1sho(index2,:) edgetoedgedata.ze1sho(index2,:);edgetoedgedata.re2sho(index2,:) edgetoedgedata.thetae2sho(index2,:) edgetoedgedata.ze2sho(index2,:)];
                    end
                    if Ndifforder >= 7
                        index1 = reftoshortlistE(edgepattern(7),edgepattern(6));
                        cylE7_r6 = [edgetoedgedata.re1sho(index1,:) edgetoedgedata.thetae1sho(index1,:) edgetoedgedata.ze1sho(index1,:);edgetoedgedata.re2sho(index1,:) edgetoedgedata.thetae2sho(index1,:) edgetoedgedata.ze2sho(index1,:)];
                        index2 = reftoshortlistE(edgepattern(6),edgepattern(7));
                        cylE6_r7 = [edgetoedgedata.re1sho(index2,:) edgetoedgedata.thetae1sho(index2,:) edgetoedgedata.ze1sho(index2,:);edgetoedgedata.re2sho(index2,:) edgetoedgedata.thetae2sho(index2,:) edgetoedgedata.ze2sho(index2,:)];
                    end
                    if Ndifforder >= 8
                        index1 = reftoshortlistE(edgepattern(8),edgepattern(7));
                        cylE8_r7 = [edgetoedgedata.re1sho(index1,:) edgetoedgedata.thetae1sho(index1,:) edgetoedgedata.ze1sho(index1,:);edgetoedgedata.re2sho(index1,:) edgetoedgedata.thetae2sho(index1,:) edgetoedgedata.ze2sho(index1,:)];
                        index2 = reftoshortlistE(edgepattern(7),edgepattern(8));
                        cylE7_r8 = [edgetoedgedata.re1sho(index2,:) edgetoedgedata.thetae1sho(index2,:) edgetoedgedata.ze1sho(index2,:);edgetoedgedata.re2sho(index2,:) edgetoedgedata.thetae2sho(index2,:) edgetoedgedata.ze2sho(index2,:)];
                    end

                    if Ndifforder == 2

                        cylE1_r2frac = cylE1_r2;
                        e1length = edgedata.edgelengthvec(edgepattern(1));
                        cylE2_r1frac = cylE2_r1;
                        e2length = edgedata.edgelengthvec(edgepattern(2));

                        [irnew,~] = EDwedge2nd(cylS,cylR,cylE2_r1frac,cylE1_r2frac,...
                            nyvec(edgepattern(:)),[0 edgedata.edgelengthvec(edgepattern(1));0 edgedata.edgelengthvec(edgepattern(2))],dzvec(edgepattern(:)),...
                            'n',pathalongplane(ii),Rstart,bc,envdata.cair,fs);

        %                 IRDIFFVEC = [IRDIFFVEC;sum(irnew)];

                    elseif Ndifforder == 3    %   if Ndifforder == 2			

                        for kk = 1:newndivvec(1)  
                            
                            BIGEDGE1stvalue = (kk-0.5)./newndivvec(1);
                            wedgeparams = [cylE2_r1;cylE1_r2;cylE3_r2;cylE2_r3];

                            [irnewpartition,~] = EDwedgeN(cylS,cylR,wedgeparams,ncylrows,...
                                nyvec(edgepattern(:)),edgedata.edgelengthvec(edgepattern(:)).',...
                                dzvec(edgepattern(:)),'n',pathalongplane(ii,:),nedgeelcombs,Rstart,bc,envdata.cair,fs,BIGEDGE1stvalue);
                            irnewpartition = real(irnewpartition);

                            if kk == 1
                                irnew = irnewpartition;    
                            else
                                lengthaddition = length(irnewpartition);
                                lengthaccum = length(irnew);
                                if lengthaddition > lengthaccum
                                    irnew = [irnew;zeros(lengthaddition-lengthaccum,1)];    
                                end
                                irnew(1:lengthaddition) = irnew(1:lengthaddition) + irnewpartition; 
                            end
                        end  
                    end
                    if Ndifforder == 4
                            wedgeparams = [cylE2_r1;cylE1_r2;cylE3_r2;cylE2_r3;cylE4_r3;cylE3_r4];                                
                    elseif Ndifforder == 5
                            wedgeparams = [cylE2_r1;cylE1_r2;cylE3_r2;cylE2_r3;cylE4_r3;cylE3_r4;cylE5_r4;cylE4_r5];
                    elseif Ndifforder == 6
                            wedgeparams = [cylE2_r1;cylE1_r2;cylE3_r2;cylE2_r3;cylE4_r3;cylE3_r4;cylE5_r4;cylE4_r5;cylE6_r5;cylE5_r6];
                    elseif Ndifforder == 7
                            wedgeparams = [cylE2_r1;cylE1_r2;cylE3_r2;cylE2_r3;cylE4_r3;cylE3_r4;cylE5_r4;cylE4_r5;cylE6_r5;cylE5_r6;cylE7_r6;cylE6_r7];
                    elseif Ndifforder == 8                               
                            wedgeparams = [cylE2_r1;cylE1_r2;cylE3_r2;cylE2_r3;cylE4_r3;cylE3_r4;cylE5_r4;cylE4_r5;cylE6_r5;cylE5_r6;cylE7_r6;cylE6_r7;cylE8_r7;cylE7_r8];
                    end

                    if Ndifforder >= 4    

                        for kk = 1:newndivvec(1)
                            if showtext >= 5
                                disp(['   ',int2str(kk),' of ',int2str(newndivvec(1))]) 
                            end
                            BIGEDGE1stvalue = (kk-0.5)./newndivvec(1);


                            [irnewpartition,~] = EDwedgeN(cylS,cylR,wedgeparams,ncylrows,...
                                nyvec(edgepattern(:)),edgedata.edgelengthvec(edgepattern(:)).',...
                                dzvec(edgepattern(:)),'n',pathalongplane(ii,:),nedgeelcombs,Rstart,bc,envdata.cair,fs,BIGEDGE1stvalue);
                            irnewpartition = real(irnewpartition);

                            if kk == 1
                                irnew = irnewpartition;
                            else
                                lengthaddition = length(irnewpartition);
                                lengthaccum = length(irnew);
                                if lengthaddition > lengthaccum
                                    irnew = [irnew;zeros(lengthaddition-lengthaccum,1)];    
                                end
                                irnew(1:lengthaddition) = irnew(1:lengthaddition) + irnewpartition; 
                            end
                        end
                    end

                    % For thin plates, we must have a boost factor!
                    % This is because there will be multiple equivalent
                    % combinations passing on the rear side of the thin plate
                    %
                    % Bug found 1 May 2017:
                    % boostfactor got erroneously set if double
                    % diffraction involved two thin planes that
                    % were not the same.                            

                    boostfactor = 1;
%                    29 May 2018: The line below caused an error for a thin
%                    plane, so it was replaced by the line below.
%                    if all( nyvec(edgepattern(:)) == 0.5 ) && pathalongplane==1                                
                   if all( nyvec(edgepattern(:)) == 0.5 ) && all(pathalongplane(ii,:)==1)                                
                        boostfactor = boostfactor*2^(Ndifforder-1);                        
                    end

                    if boostfactor ~= 1
                        irnew = irnew*boostfactor;    
                    end

                    nnew = length(irnew);

                    % The case with size(irhod,2) < cellcounter
                    % could happen if an S/R combo larger than numbers 1/1
                    % has a higher difforder than S/R = 1/1.
                    
                    if firstcomponentdone == 0 || size(irhod,2) < cellcounter
                       if Sdata.doaddsources == 0
                           irhod{cellcounter} = zeros(nnew,nreceivers,nsources);
                           irhod{cellcounter}(:,irec,isou) = irnew*Sdata.sourceamplitudes(isou);
                       else
                           irhod{cellcounter} = zeros(nnew,nreceivers,1);
                           irhod{cellcounter}(:,irec,1) = irnew*Sdata.sourceamplitudes(isou);                   
                       end
                       firstcomponentdone = 1;
                    else
                        nold = size(irhod{cellcounter},1);
                        if nnew > nold
                            if Sdata.doaddsources == 0
                               irhod{cellcounter} = [irhod{cellcounter};zeros(nnew-nold,nreceivers,nsources)]; 
                            else
                               irhod{cellcounter} = [irhod{cellcounter};zeros(nnew-nold,nreceivers,1)];                         
                            end
                        end
                        if Sdata.doaddsources == 0
                            irhod{cellcounter}(1:nnew,irec,isou) = irhod{cellcounter}(1:nnew,irec,isou) + irnew*Sdata.sourceamplitudes(isou);                
                        else
                            irhod{cellcounter}(1:nnew,irec,1) = irhod{cellcounter}(1:nnew,irec,1) + irnew*Sdata.sourceamplitudes(isou);                                    
                        end
                    end

                end 
            else
                disp(['   No paths for order ',int2str(Ndifforder)])
            end  

        end  

    end
    
end

if controlparameters.savealldifforders == 1
    nirlength = 0;
   for ii = 2:difforder
       nnew = size(irhod{ii},1);
       if nnew > nirlength
           nirlength = nnew;
       end
   end
   
   for ii = 2:difforder
       nnew = size(irhod{ii},1);
       if nnew < nirlength
           if Sdata.doaddsources == 1
               irhod{ii} = [irhod{ii};zeros(nirlength-nnew,nreceivers)];
           else
               irhod{ii} = [irhod{ii};zeros(nirlength-nnew,nreceivers,nsources)];               
           end
       end
   end
      
end

elapsedtimemakeirhod = etime(clock,t00);

save(desiredname,'irhod','EDinputdatahash','elapsedtimemakeirhod','stageversion','EDversionnumber');
