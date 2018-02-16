function pathstruct = EDfindpathsISESx(planedata,edgedata,lengthNspecmatrix,lengthNdiffmatrix,...
    singlediffcol,startindicessinglediff,endindicessinglediff,...
    ndecimaldivider,PointertoIRcombs,IRoriginsfrom,S,R,...
    directsound,specorder,difforder,nedgesubs,visplanesfromS,visplanesfromR,...
    vispartedgesfromS,vispartedgesfromR,showtext)
% EDfindpathsISESx - Finds all possible paths that include direct sound, specular, diffraction.
% Finds all possible paths of types direct sound, specular, diffraction,
% combinations of specular and diffraction. This is done by going through
% the potentially visible set in the ISEStree. The data of this tree is
% passed as global variables.
%
% Input parameters:
%   planedata,edgedata  structs
%   lengthNspecmatrix,lengthNdiffmatrix,...
%   singlediffcol,startindicessinglediff,endindicessinglediff,...
%   ndecimaldivider,PointertoIRcombs,IRoriginsfrom     Vectors from the
%                       ISEStree
% 	S, R,               Coordinates and counter number of source and receiver, taken
%                       from the setup file.
% 	directsound         0 or 1
%   specorder           0 or higher integer. The desired max. refl. order
%   difforder           0 or higher integer, <= specorder. The desired max.
%                       diffr. order
% 	nedgesubs           The number of divisions of each edge, for
%                       a simplified edge visibility check.
% 	visplanesfromS, visplanesfromR, vispartedgesfromS, vispartedgesfromR
%                       Taken directly from the srdatafiles.
%   showtext
% Global paramaters:
%   IVNDIFFMATRIX       A list from the ISEStreefile
% 
% These global parameters from the ISEStreefile are accessed by subroutines
%   ISCOORDS POTENTIALISES ORIGINSFROM ISESVISIBILITY 
%   IVNSPECMATRIX REFLORDER
%
% Output parameters:
%   pathstruct              struct, with these fields:
%       S                       The coordinates of the source.
%       R                       The coordinates of the receiver.
%       Sinsideplanenumber      The plane number that the source is directly at. If the
%                           source is not placed directly at a plane, then
%                           souinsideplanenumber = [];
%       Rinsideplanenumber      The plane number that the receiver is directly at. If the
%                           receiver is not placed directly at a plane, then
%                           recinsideplanenumber = [];
%       mainlistguide           A matrix, [nposscombs,3], which for each row gives
%                           1. the number of rows in "reflpaths" and
%                           "pathtypevec" that have one type of combination.
%                           2. the first row and 3. the last row
%       mainlistguidepattern    A matrix, [nposscombs,specorder] which for each
% 	                        row gives the description of each possible combination,
%                           using 'f', 's' and 'd'.
%       directsoundrow          The value 0 or 1 indicating whether the first row in
%                           mainlistguide & mainlistguidepattern & reflpaths & pathtypevec 
%                           contains the direct sound or not.
%       allspecrows             A list, [1,2], containing the first and last row numbers
%                           of mainlistguide & mainlistguidepattern that contain specular
%                           reflections. If there are no specular reflections, allspecrows
%                           contains [0 0].
%       firstdiffrow            A number indicating which is the first row of mainlistguide &
%                           mainlistguidepattern that contain a diffraction component. If
%                           there are no diffraction components, firstdiffrow = 0.
%       pathtypevec             A matrix, [ncombs,specorder], describing the
% 	                        type of reflection, using 'f' for the direct
% 	                        sound, 's' for specular reflection and 'd' for
% 	                        diffraction.
%       reflpaths               A matrix, [ncombs,specorder], giving the plane 
%                           and edge numbers involved in each refl. comb.
%       specextradata           A sparse matrix, [ncombs,(specorder+1)*3],
% 	                        containing the coordinates of the image source
% 	                        (col 1-3), the image receiver (col 4-6)
%                           and the coordinates of all specular hit points
%                           (only for purely specular combinations).
%       edgeextradata           A sparse matrix, [ncombs,(specorder*2)], with
% 	                        the visibility of the involved edge, denoted
%                           by two values between 0 and 1.
%
% Uses the functions EDdirectsound EDspeculISES EDdiffISESx
%                    ESIE2diff2ISES ESIE2diffNISES
%
% ----------------------------------------------------------------------------------------------
% Peter Svensson (peter.svensson@ntnu.no) 28 Nov. 2017
%
% [pathtypevec,reflpaths,specextradata,edgeextradata,S,R,mainlistguide,mainlistguidepattern,...
%     directsoundrow,allspecrows,firstdiffrow,Sinsideplanenumber,Rinsideplanenumber] = ...
%     EDfindpathsISESx(planedata,edgedata,...
%     lengthNspecmatrix,lengthNdiffmatrix,singlediffcol,startindicessinglediff,...
%     endindicessinglediff,ndecimaldivider,PointertoIRcombs,IRoriginsfrom,S,R,isou,irec,...
%     directsound,specorder,difforder,nedgesubs,visplanesfromS,visplanesfromR,...
%     vispartedgesfromS,vispartedgesfromR,showtext);

% 6 Apr 2017    Reintroduced the possibility to save pathfile
% 1 Nov. 2017   Re-removed the possibility to save pathfile. Saved
%               externally to this function.
% 28 Nov. 2017  Copied to EDtoolbox. Introduced the non-global input parameter
%               showtext. Packed the output data into a struct.

global IVNDIFFMATRIX

pathtypevec = [];
reflpaths = [];
specextradata = [];
edgeextradata = [];
mainlistguide = [];
mainlistguidepattern = [];

%-------------------------------------------------------
%		#############################
%		##     DIRECT SOUND        ##
%		#############################

if directsound == 1
	if showtext >= 2
		disp('      Direct sound')
	end
    dirsoundok = EDdirectsound(planedata.corners,planedata.planecorners,planedata.planeeqs,planedata.planeeqs(:,1:3),...
        S,R,visplanesfromS,visplanesfromR,planedata.canplaneobstruct,planedata.minvals,planedata.maxvals,planedata.ncornersperplanevec,showtext);
    if dirsoundok == 1
        pathtypevec = [pathtypevec;['f',zeros(1,specorder-1)]];
        reflpaths = [reflpaths;zeros(1,specorder)];
        if specorder > 0
            specextradata = [specextradata;[S R zeros(1,(specorder-1)*3)]];
        else
            specextradata = [specextradata;[S R]];   
        end
        edgeextradata = [edgeextradata;zeros(1,specorder*2)];    
        mainlistguide = [1 1 1];   
        mainlistguidepattern = ['f' ' '*ones(1,specorder-1)];
        directsoundrow = 1;
    else
        directsoundrow = 0;    
    end
else
    directsoundrow = 0;
end

%-------------------------------------------------------
%		#############################
%		##        SPECULAR         ##
%		#############################

if specorder >= 1
	if showtext >= 2
		disp('      SPECULAR reflections')
    end
    
	[validISlist,validIScoords,allreflpoints,listguide,listofreflorder] = EDspeculISES(planedata.corners,planedata.planecorners,planedata.planeeqs,planedata.planeeqs(:,1:3),...
        S,R,lengthNspecmatrix,specorder,visplanesfromR,planedata.planeisthin,planedata.canplaneobstruct,planedata.minvals,planedata.maxvals,planedata.ncornersperplanevec,...
        planedata.planeseesplane,planedata.rearsideplane,planedata.modeltype,showtext);
    nvalidreflorders = size(listguide,1);
    for ii = 1:nvalidreflorders
        ncombs = listguide(ii,1);
        norder = listofreflorder(ii);
        
        if ncombs > 0
            pathtypevec = [pathtypevec;['s'*ones(ncombs,norder) zeros(ncombs,specorder-norder)]];
            reflpaths = [reflpaths;validISlist(listguide(ii,2):listguide(ii,3),:)];
            if size(pathtypevec,2) > size(reflpaths,2)
                reflpaths = [reflpaths zeros(size(reflpaths,1),size(pathtypevec,2)-size(reflpaths,2))];    
            end
            specextradata = [specextradata;[validIScoords(listguide(ii,2):listguide(ii,3),1:3) allreflpoints(listguide(ii,2):listguide(ii,3),:)]];
            if size(specextradata,2) < (size(pathtypevec,2)+1)*3
                specextradata = [specextradata zeros(size(specextradata,1),(size(pathtypevec,2)+1)*3-size(specextradata,2))];
            end
            edgeextradata = [edgeextradata;zeros(ncombs,specorder*2)];    
        end            

        mainlistguidepattern = [mainlistguidepattern;['s'*ones(1,norder) ' '*ones(1,specorder-norder)]];
    
    end

    nprevious = size(mainlistguide,1);
    if nprevious > 0
        listguide(:,2:3) = listguide(:,2:3)+mainlistguide(nprevious,3);
    end
    mainlistguide = [mainlistguide;listguide];
    if nvalidreflorders == 0
        allspecrows = [0 0];
    else        
        allspecrows = [min([1 nvalidreflorders]) nvalidreflorders]+directsoundrow;
    end
    firstdiffrow = nvalidreflorders + directsoundrow + 1;
   
else
    allspecrows = [0 0];        
    firstdiffrow = directsoundrow + 1;
end

%-------------------------------------------------------
%		####################################
%		##      SINGLE DIFFRACTION        ##
%		####################################

if difforder >=1 && ~isempty(lengthNdiffmatrix)
	if showtext >= 2
		disp('      SINGLE DIFFRACTION combined with any order specular')
	end
    [edgedifflist,startandendpoints,prespeclist,postspeclist,validISEDcoords,validEDIRcoords,listguide,listoforders,...
        bigedgeweightlist] = EDdiffISESx(planedata,edgedata,S,R,...
        IVNDIFFMATRIX(1:lengthNdiffmatrix(1),1),singlediffcol,startindicessinglediff,endindicessinglediff,...
        specorder,visplanesfromR,vispartedgesfromS,vispartedgesfromR,nedgesubs,ndecimaldivider,PointertoIRcombs,IRoriginsfrom,showtext);    
    nvalidcombs = size(listguide,1);
    nprespecs = listoforders(:,2) - 1;
    npostspecs = listoforders(:,1) - nprespecs - 1;
    
    [nrows,nsubsegmentcols] = size(startandendpoints);
% % %     nsubsegments = (startandendpoints(:,1:2:nsubsegmentcols)+startandendpoints(:,2:2:nsubsegmentcols))~=0;
    nsubsegments = uint8( (startandendpoints(:,2:2:nsubsegmentcols))~=0 );

    if nrows > 1 && nsubsegmentcols > 2
        nsubsegments = sum(nsubsegments.').';
    end
    nmaxsubsegments = max(nsubsegments);
    
    expandedlistguide = listguide;
    for ii = 1:nvalidcombs
        ncombs = listguide(ii,1);
        iv = [listguide(ii,2):listguide(ii,3)];

        pathtypevec   = [pathtypevec;  ['s'*ones(ncombs,nprespecs(ii))       'd'*ones(ncombs,1)         's'*ones(ncombs,npostspecs(ii))       zeros(ncombs,specorder-listoforders(ii)) ]];
        reflpaths     = [reflpaths;    [prespeclist(iv,1:nprespecs(ii))   edgedifflist(iv)         postspeclist(iv,1:npostspecs(ii))  zeros(ncombs,specorder-listoforders(ii)) ]];  
        specextradata = [specextradata;[  validISEDcoords(iv,1:3)         validEDIRcoords(iv,1:3)    zeros(ncombs,(specorder-1)*3)            ]];
        edgeextradata = [edgeextradata;[startandendpoints(iv,1:2)           zeros(ncombs,(specorder-1)*2)                         ]];    
        mainlistguidepattern = [mainlistguidepattern;['s'*ones(1,nprespecs(ii))       'd'         's'*ones(1,npostspecs(ii))       zeros(1,specorder-listoforders(ii))]];
        if nmaxsubsegments > 2
        for jj = 2:nmaxsubsegments
            ivsub = find(nsubsegments(iv)>=jj);
            ncombssubset = length(ivsub);
            if ncombssubset > 0
                pathtypevec   = [pathtypevec;  ['s'*ones(ncombssubset,nprespecs(ii))       'd'*ones(ncombssubset,1)         's'*ones(ncombssubset,npostspecs(ii))       zeros(ncombssubset,specorder-listoforders(ii)) ]];
                reflpaths     = [reflpaths;    [prespeclist(iv(ivsub),1:nprespecs(ii))   edgedifflist(iv(ivsub))         postspeclist(iv(ivsub),1:npostspecs(ii))  zeros(ncombssubset,specorder-listoforders(ii)) ]];  
                specextradata = [specextradata;[  validISEDcoords(iv(ivsub),1:3)         validEDIRcoords(iv(ivsub),1:3)    zeros(ncombssubset,(specorder-1)*3)            ]];
                edgeextradata = [edgeextradata;[startandendpoints(iv(ivsub),(jj-1)*2+1:(jj-1)*2+2)           zeros(ncombssubset,(specorder-1)*2)                         ]];
                expandedlistguide(ii,1) = expandedlistguide(ii,1) + ncombssubset;
                expandedlistguide(ii,3) = expandedlistguide(ii,3) + ncombssubset;
                for kk = ii+1:nvalidcombs
                    expandedlistguide(kk,2:3) = expandedlistguide(kk,2:3) + ncombssubset; 
                end
            end    
        end
        end
    end
    n1 = size(mainlistguide,1);
    if n1 > 0
        listoffset = mainlistguide(n1,3);
    else
        listoffset = 0;    
    end
    if size(expandedlistguide,1) > 0
        mainlistguide = [mainlistguide;[expandedlistguide(:,1) expandedlistguide(:,2:3)+listoffset]];
    end
end

%-------------------------------------------------------
%		################################################
%		##      DOUBLE AND HIGHER DIFFRACTION         ##
%		################################################

for kk = 2:min([specorder difforder])

    if ~isempty(lengthNdiffmatrix)
        
    if showtext >= 2
		disp(['   DIFFRACTION, ORDER ',int2str(kk),' combined with any order specular'])
	end

    nmaxsubsegments = 0;
    
    % We will use the spec-edge-spec combs that were found valid for
    % first-order diffraction, because if higher-order diffraction combs
    % end with the same edge-spec-spec-... comb we know the visibility
    % already.
    % First, pick out only the combs that have postspecs. Second, keep
    % only a unique set because there will be many repeated combinations.
    
   
    if kk == 2
        iv = find(sum(postspeclist.')>0 );
        edgedifflist = edgedifflist(iv,:);
        postspeclist = postspeclist(iv,:);
        bigedgeweightlist = bigedgeweightlist(iv,:);
        validEDIRcoords = validEDIRcoords(iv,:);

        patternthatisOK = [edgedifflist postspeclist];
        [~,iv,~] = unique(patternthatisOK,'rows');
        edgedifflistin = edgedifflist(iv,:);
        postspeclistin = postspeclist(iv,:);
        bigedgeweightlistin = bigedgeweightlist(iv,:);
        validEDIRcoordsin = validEDIRcoords(iv,:);

    end
    
    if kk == 2
        maxrownumber = max(lengthNdiffmatrix(1:2));
        [edgedifflist,startandendpoints,prespeclist,midspeclist,postspeclist,validISEDcoords,validEDIRcoords,listguide,listofallspecs] = ESIE2diff2ISES(planedata,edgedata,S,R,...
            IVNDIFFMATRIX(1:maxrownumber,1:2),lengthNdiffmatrix(1:2),...
            specorder,visplanesfromR,vispartedgesfromS,vispartedgesfromR,nedgesubs,ndecimaldivider,edgedifflistin,postspeclistin,...
            bigedgeweightlistin,validEDIRcoordsin);

        [nrows,nsubsegmentcols] = size(startandendpoints);
        nsubsegments = uint8((startandendpoints(:,2:4:nsubsegmentcols))~=0);
        if nrows > 1 && nsubsegmentcols > 4
            nsubsegments = sum(nsubsegments.').';
        end
        nmaxsubsegments = max(nsubsegments);

    elseif kk >= 3
        maxrownumber = max(lengthNdiffmatrix(1:kk));
        [edgedifflist,startandendpoints,prespeclist,postspeclist,validISEDcoords,validEDIRcoords,listguide,listoforders] = ESIE2diffNISES(planedata,edgedata,S,R,...
            IVNDIFFMATRIX(1:maxrownumber,1:kk),lengthNdiffmatrix(1:kk),kk,...
            specorder,visplanesfromR,vispartedgesfromS,vispartedgesfromR,nedgesubs,ndecimaldivider,edgedifflistin,postspeclistin,...
            bigedgeweightlistin,validEDIRcoordsin);
    end

    [nvalidcombs,~] = size(listguide); 
    if kk >= 3
        nprespecs = listoforders(:,2) - 1;
        npostspecs = listoforders(:,1) - nprespecs - kk;
    else
        nprespecs = listofallspecs(:,1);    
        nmidspecs = listofallspecs(:,2);    
        npostspecs = listofallspecs(:,3);    
        listoforders = nprespecs+nmidspecs+npostspecs+2;
    end

        
    for ii = 1:nvalidcombs
        ncombs = listguide(ii,1);
        i1 = listguide(ii,2);     i2 = listguide(ii,3);
        if kk == 2
            pathtypevec   = [pathtypevec;  ['s'*ones(ncombs,nprespecs(ii))  'd'*ones(ncombs,1)   's'*ones(ncombs,nmidspecs(ii))  'd'*ones(ncombs,1)    's'*ones(ncombs,npostspecs(ii))       zeros(ncombs,specorder-listoforders(ii)) ]];
        else
            pathtypevec   = [pathtypevec;  ['s'*ones(ncombs,nprespecs(ii))  'd'*ones(ncombs,kk)   's'*ones(ncombs,npostspecs(ii))       zeros(ncombs,specorder-listoforders(ii)) ]];            
        end
        
        if kk == 2        
            reflpaths     = [reflpaths;    [prespeclist(i1:i2,1:nprespecs(ii))   edgedifflist(i1:i2,1)  midspeclist(i1:i2,1:nmidspecs(ii)) edgedifflist(i1:i2,2)         postspeclist(i1:i2,1:npostspecs(ii))  zeros(ncombs,specorder-listoforders(ii)) ]];  
        else
            reflpaths     = [reflpaths;    [prespeclist(i1:i2,1:nprespecs(ii))   edgedifflist(i1:i2,1:kk)         postspeclist(i1:i2,1:npostspecs(ii))  zeros(ncombs,specorder-listoforders(ii)) ]];              
        end
        
        specextradata = [specextradata;[  validISEDcoords(i1:i2,1:3)         validEDIRcoords(i1:i2,1:3)    zeros(ncombs,(specorder-1)*3)            ]];
        edgeextradata = [edgeextradata;[startandendpoints(i1:i2,1:kk*2)           zeros(ncombs,(specorder-kk)*2)                         ]];    
        if kk == 2
            mainlistguidepattern = [mainlistguidepattern;['s'*ones(1,nprespecs(ii))   'd' 's'*ones(1,nmidspecs(ii)) 'd'         's'*ones(1,npostspecs(ii))       zeros(1,specorder-listoforders(ii))]];
        else
            mainlistguidepattern = [mainlistguidepattern;['s'*ones(1,nprespecs(ii))   'd'*ones(1,kk)  's'*ones(1,npostspecs(ii))       zeros(1,specorder-listoforders(ii))]];            
        end
        
        for jj = 2:nmaxsubsegments
            ivsub = find(nsubsegments(iv)>=jj);
            ncombssubset = length(ivsub);
            if ncombssubset > 0
                pathtypevec   = [pathtypevec;  ['s'*ones(ncombssubset,nprespecs(ii))       'd'*ones(ncombssubset,1) 's'*ones(ncombssubset,nmidspecs(ii)) 'd'*ones(ncombssubset,1)        's'*ones(ncombssubset,npostspecs(ii))       zeros(ncombssubset,specorder-listoforders(ii)) ]];
                reflpaths     = [reflpaths;    [prespeclist(iv(ivsub),1:nprespecs(ii))   edgedifflist(iv(ivsub),1)  midspeclist(iv(ivsub),1:nmidspecs(ii))  edgedifflist(iv(ivsub),2)         postspeclist(iv(ivsub),1:npostspecs(ii))  zeros(ncombssubset,specorder-listoforders(ii)) ]];  
                specextradata = [specextradata;[  validISEDcoords(iv(ivsub),1:3)         validEDIRcoords(iv(ivsub),1:3)    zeros(ncombssubset,(specorder-1)*3)            ]];
                edgeextradata = [edgeextradata;[startandendpoints(iv(ivsub),(jj-1)*4+1:(jj-1)*4+4)           zeros(ncombssubset,(specorder-1)*4)                         ]];
                expandedlistguide(ii,1) = expandedlistguide(ii,1) + ncombssubset;
                expandedlistguide(ii,3) = expandedlistguide(ii,3) + ncombssubset;
                for kk = ii+1:nvalidcombs
                    expandedlistguide(kk,2:3) = expandedlistguide(kk,2:3) + ncombssubset; 
                end
            end    
         end
         
     end

    n1 = size(mainlistguide,1);
    if n1 > 0
        listoffset = mainlistguide(n1,3);
    else
        listoffset = 0;    
    end
    mainlistguide = [mainlistguide;[listguide(:,1) listguide(:,2:3)+listoffset]];

    end

end    

%-------------------------------------------------------

n1 = size(mainlistguide,1);
if firstdiffrow > n1
    firstdiffrow = 0;    
end

%-------------------------------------------------------
% Is the source or receiver directly at a plane?

Sinsideplanenumber = find(visplanesfromS==4);

Rinsideplanenumber = find(visplanesfromR==4);

%-------------------------------------------------------
% Get rid of empty columns

if ~isempty(pathtypevec)
    checksum = sum(pathtypevec);
    ncols = length(checksum);
    while checksum(ncols) == 0
        pathtypevec(:,ncols) = [];
        checksum = sum(pathtypevec);
        ncols = length(checksum);
    end
end
    
%-------------------------------------------------------
% Compact the output data

specextradata = sparse(specextradata);
edgeextradata = sparse(edgeextradata);
pathtypevec = uint8(pathtypevec);
ncombs = size(reflpaths,1);
if ncombs+1 < 256
    mainlistguide = uint8(mainlistguide);
elseif ncombs+1 < 65536    
    mainlistguide = uint16(mainlistguide);
else
    mainlistguide = uint32(mainlistguide);
end

pathstruct = struct('S',S);
pathstruct.R = R;
pathstruct.Sinsideplanenumber = Sinsideplanenumber;
pathstruct.Rinsideplanenumber = Rinsideplanenumber;
pathstruct.mainlistguide = mainlistguide;
pathstruct.mainlistguidepattern = mainlistguidepattern;
pathstruct.directsoundrow = directsoundrow;
pathstruct.allspecrows = allspecrows;
pathstruct.firstdiffrow = firstdiffrow;
pathstruct.pathtypevec = pathtypevec;
pathstruct.reflpaths = reflpaths;
pathstruct.specextradata = specextradata;
pathstruct.edgeextradata = edgeextradata;

