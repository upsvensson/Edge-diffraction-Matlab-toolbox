function [edgedifflist,startandendpoints,prespeclist,postspeclist,validIScoords,validIRcoords,listguide,listoforders,...
        bigedgeweightlist] = EDdiffISESx(planedata,edgedata,S,R,...
        ivsinglediff,singlediffcol,startindicessinglediff,endindicessinglediff,...
        specorder,visplanesfromR,vispartedgesfromS,vispartedgesfromR,nedgesubs,ndecimaldivider,...
        PointertoIRcombs,IRoriginsfrom,showtext)
% EDdiffISESx - Gives list of paths that includes a 1st-order diff. combination.
% Gives the list of possible first-order diffraction paths, possibly with specular
% reflections before and after.
%
% Input parameters:
%   planedata,edgedata  structs
%   S,R
%   ivsinglediff,singlediffcol,startindicessinglediff,endindicessinglediff
%                   from the ISEStreefile
%   specorder       The max. wanted reflection order
%   visplanesfromR,vispartedgesfromS,vispartedgesfromR
%   nedgesubs       The number of edge divisions for simplified edge
%                   visibility check
%   ndecimaldivider,PointertoIRcombs,IRoriginsfrom
%                   from the ISEStreefile
%   showtext
%
% Global parameters:
%   POTENTIALISES ISCOORDS ORIGINSFROM REFLORDER ISESVISIBILITY
%                       From the ISEStreefile
%
% Output parameters:
%   edgedifflist        List [ncombs,1] of the edge number involved in each
%                       spec-diff-spec combination.
%   startandendpoints   Matrix [ncombs,2] of the relative start and end
%                       points of each edge. The values, [0,1], indicate
%                       which part of the edge that is visible.
%   prespeclist         Matrix [ncombs,specorder-1] of the specular
%                       reflections that precede every diffraction.
%   postspeclist        Matrix [ncombs,specorder-1] of the specular
%                       reflections that follow every diffraction.
%   validIScoords       Matrix [ncombs,3] of the image source for each
%                       multiple-spec that precedes the diffraction. If
%                       there is no spec refl before the diffraction, the
%                       value [0 0 0] is given.
%   validIRcoords       Matrix [ncombs,3] of the image receiver for each
%                       multiple-spec that follows the diffraction. If
%                       there is no spec refl after the diffraction, the
%                       value [0 0 0] is given.
%   listguide           Matrix [nuniquecombs,3] which for each row gives
%                       1. The number of examples in edgefdifflist etc that
%                          are the same type of spec-diff-spec comb.
%                       2. The first row number and 3. The last row number.
%   listoforders        Matrix [nuniquecombs,2] which for each row gives
%                       1. The reflection order for the spec-diff-spec comb
%                          in the same row in listguide.
%                       2. The order of the diffraction in the
%                          spec-diff-spec comb.
%   bigedgeweightlist   List [ncombs,1] of the visibility of each edge,
%                       expressed as a number 0 to 2^nedgesubs-1.
%
% Uses functions  EDfindis EDgetedgepoints EDchkISvisible EDcheckobstrpaths
%
% Peter Svensson (peter.svensson@ntnu.no) 16 March 2021
%
% [edgedifflist,startandendpoints,prespeclist,postspeclist,validIScoords,validIRcoords,listguide,listoforders,...
%   bigedgeweightlist] = EDdiffISESx(planedata,edgedata,S,R,...
%         ivsinglediff,singlediffcol,startindicessinglediff,endindicessinglediff,...
%         specorder,visplanesfromR,vispartedgesfromS,vispartedgesfromR,nedgesubs,ndecimaldivider,...
%         PointertoIRcombs,IRoriginsfrom)

% 18 Nov. 2006 Functioning version
% 28 Nov. 2017 Copied to EDtoolbox. Introduced the showtext non-global
% input parameter
% 16 Mar 2021 Adapted to changes in EDchkISvisible.

global POTENTIALISES ISCOORDS ORIGINSFROM REFLORDER ISESVISIBILITY

%  eval(['load ',eddatafile])
clear edgeseesplane cornerinfrontofplane

nedges = size(edgedata.planesatedge,1);
nplanes = size(planedata.planecorners,1);
planenvecs = planedata.planeeqs(:,1:3);

edgedifflist = [];
prespeclist = [];
postspeclist = [];
startandendpoints = [];
bigedgeweightlist = [];
validIScoords = [];
validIRcoords = [];
%%%%%allreflpoints = [];

n2 = size(POTENTIALISES,2);
if n2 < specorder
    specorder = n2;    
end

maxvisibilityvalue = 2^nedgesubs-1;
zerosvec1 = zeros(1,max([specorder-1 1]));
% zerosvec2 = zeros(1,3);
% listguide = zeros(specorder*2-1,3);
listoforders = zeros(specorder*2-1,2);

obstructtestneeded = (sum(planedata.canplaneobstruct)~=0);
% onesvec = ones(1,nedgesubs);
% onesvec3 = ones(1,3);

%   ###########################################
%   #                                         #
%   #         S - spec-spec- edge - R cases   #
%   #                                         #
%   #         Prespec cases                   #
%   #                                         #
%   ###########################################
%
% Possible edges for S-spec-spec-E-R are seen (at least partly) by the receiver.
%
% The visibility doesn't need to be checked since this was done in the
% ISEStree.

% The vector masterivlist will always refer to the original data vector,
% i.e. POTENTIALISES, ORIGINSFROM, REFLORDER etc
%
% First we pick out those indices where there was a single diffraction, but
% skip those with only diffraction (because we dealt with them already).
% Also, select only those where the diffraction is the last in the sequence
% of reflections.

for ii = 1:specorder      % NB!!! ii = 1 corresponds to zero spec. refl before the diff.
                            %       ii = 2 corresponds to one spec refl.
                            %       before the diff etc

    if showtext >= 3
        if ii > 1
            disp(['      Checking for ',int2str(ii-1),' spec refl before the edge diff'])    
        else
            disp('      Checking for 0 spec refl before the edge diff')    
        end
    end

    iv = uint32(startindicessinglediff(ii):endindicessinglediff(ii));
        
    if ~isempty(iv)
                
        ivkeep = singlediffcol(iv)==ii;
        iv = iv(ivkeep);
        masterivlist = ivsinglediff(iv);
        possibleedges = double(POTENTIALISES(masterivlist,ii)) - nplanes;

        % Keep only combinations for which the receiver can see the edge
            
        ivnotvisiblefromr = find(vispartedgesfromR(possibleedges)==0);
        if ~isempty(ivnotvisiblefromr)
            masterivlist(ivnotvisiblefromr) = [];
            possibleedges(ivnotvisiblefromr) = [];
        end        
        % Pick out the pre-specs

        nposs = length(masterivlist);

        if nposs > 0
        
            if ii > 1
                possiblecombs = POTENTIALISES(masterivlist,1:ii-1);
%                 reftoIScoords = ORIGINSFROM(masterivlist);
                edgeweightlist = bitand((ISESVISIBILITY(masterivlist)),(vispartedgesfromR(possibleedges)));
                prespeclist = [prespeclist;[possiblecombs zeros(nposs,specorder-ii)]];
                % NB! It is correct below that the indices for the IScoords should be
                % ORIGINSFROM(masterivlist), rather than masterivlist.
                % The combinations in POTENTIALISES(masterivlist,:) all have
                % spec-spec-...-diff combinations and then
                % ISCOORDS(masterivlist,:) are zeros since a comb. that
                % ends with a diff has no image source. 
                 validIScoords = [validIScoords;ISCOORDS(ORIGINSFROM(masterivlist),:)];
            else            
                edgeweightlist = bitand((vispartedgesfromS(possibleedges)),(vispartedgesfromR(possibleedges)));
                prespeclist = [prespeclist;[zeros(nposs,max([specorder-1 1]))]];       
                % For the case of no spec refl before the diffraction, we
                % let the IS get the S coordinates.
                validIScoords = [validIScoords;S(ones(nposs,1),:)];
            end
            
            if showtext >= 3
         		disp(['         ',int2str(nposs),' IS+edges valid'])
            end
		
            edgedifflist = [edgedifflist;possibleedges];
            postspeclist = [postspeclist;zerosvec1(ones(nposs,1),:)];
            bigedgeweightlist = [bigedgeweightlist;edgeweightlist];   
            validIRcoords = [validIRcoords;R(ones(nposs,1),:)];
        
        end
        
    end
end

iv = uint32(find(bigedgeweightlist==0));
if ~isempty(iv)
    edgedifflist(iv) = [];
    bigedgeweightlist(iv) = [];
    prespeclist(iv,:) = [];
    postspeclist(iv,:) = [];
    validIScoords(iv,:) = [];
    validIRcoords(iv,:) = [];
end

%   #######################################################################
%   #######################################################################
%   #######################################################################
%   #######################################################################
%   #######################################################################

%   #############################################
%   #                                           #
%   #         Build an IR tree                 #
%   #                                           #
%   #############################################
%
% For all cases with a specular reflection after an edge diffraction,
% we build an IR tree, analogous to the IS tree

if specorder > 1

    if showtext >= 3
            disp('      For the edge+spec combs, an IR tree is built')    
    end

    % Select all combinations with a single diffraction, anywhere except in
    % the last column, and store the useful data for these (which column
    % is the diffraction, what reflection order, how many specular
    % reflections?)
   
    ivselect = uint32(find( singlediffcol~=REFLORDER(ivsinglediff) ));
    masterivlistorig = ivsinglediff(ivselect);
    singlediffcol_select = singlediffcol(ivselect);
    reflorder_select = REFLORDER(ivsinglediff(ivselect));
    clear ivsinglediff
    nspecular_select = uint32(double(reflorder_select) - double(singlediffcol_select));
    nspecular_select_not_empty = 1;
    
    % Now we remove all those with a last reflection plane which can not be
    % seen by the receiver
    
    if nplanes + nedges < 65536
        lastreflplane = uint16(zeros(length(ivselect),1));
    else
        lastreflplane = uint32(zeros(length(ivselect),1));        
    end
    clear ivselect
    for ii = 2:specorder
        iv = find(reflorder_select==ii);
        lastreflplane(iv) = POTENTIALISES(masterivlistorig(iv),ii);
    end
    ivremove = uint32(find(visplanesfromR(lastreflplane)~=2 & visplanesfromR(lastreflplane)~=4 & visplanesfromR(lastreflplane)~=5));
    if ~isempty(ivremove)
        masterivlistorig(ivremove) = [];
%         singlediffcol_select(ivremove) = [];
        reflorder_select(ivremove) = [];
        nspecular_select(ivremove) = [];
%         lastreflplane(ivremove) = [];
        clear ivremove
    end  
    	
    % Start by calculating all the first-order IR coordinates for
    % all the planes that are visible from the receiver
	IRcoordsstart = zeros(nplanes,3);
	iv = uint32(find(visplanesfromR==2 | visplanesfromR==4 | visplanesfromR==5));
	IRcoordsstart(iv,:) = EDfindis(R,iv,planedata.planeeqs);
		
	bigIRcoords = (zeros(size(ISCOORDS)));
% 	IRreflist = zeros(size(iv));
	
	% IMPROVE: Not necessary to calculate the IR coordinates for all
	% combinations since many are the same!

    % Now we make a temporary copy of POTENTIALISES which is displaced to
    % the right since this will make it easier to align the postspec
    % combos.
    
    PotentialISESshift = POTENTIALISES;
    for ii = 2:specorder-1
        iv = uint32(find(reflorder_select==ii));    
        PotentialISESshift(masterivlistorig(iv),:) = [zeros(length(iv),specorder-ii) PotentialISESshift(masterivlistorig(iv),1:ii)];
    end
    
    % Go through ii, which is the number of specular reflections after
    % the diffraction.

    if ~isempty(nspecular_select)
		for ii = 1:specorder-1
			if showtext >= 3
                    disp(['         order ',int2str(ii)])    
			end
		
            % We save the index vectors from the different orders so they can
            % be reused in the next stage.
			
             iv = uint32(find(nspecular_select == ii));
             listselect = masterivlistorig(iv);
%              nlist = length(listselect);
             eval(['iv',int2str(ii),' = iv;']) 
                 
            if ii == 1
                bigIRcoords(listselect,:) = IRcoordsstart( PotentialISESshift(listselect,specorder),:);
        
            elseif ii == 2
                bigIRcoords(listselect,:) = EDfindis(IRcoordsstart(PotentialISESshift(listselect,specorder),:),...
                                                               PotentialISESshift(listselect,specorder-1),planedata.planeeqs);
	
             elseif ii >= 3
                ivcombo = double(PotentialISESshift(listselect,specorder));
                for jj = 1:ii-2
                    ivcombo = ivcombo + double(PotentialISESshift(listselect,specorder-jj))*ndecimaldivider^(jj);              
                end
                IRreflist = PointertoIRcombs(ivcombo);
                ivsimple    = find(IRreflist~=0);
                ivnotsimple = find(IRreflist==0);
                newIRcoordssimple    = EDfindis(bigIRcoords(IRreflist(ivsimple),:),...
                                                  PotentialISESshift(listselect(ivsimple),specorder-(ii-1)),planedata.planeeqs);
                newIRcoordsnotsimple = EDfindis(IRcoordsstart(PotentialISESshift(listselect(ivnotsimple),specorder),:),...
                                                  PotentialISESshift(listselect(ivnotsimple),specorder-1),planedata.planeeqs);
                for jj = 2:ii-1
                    newIRcoordsnotsimple = EDfindis(newIRcoordsnotsimple,...
                                                  PotentialISESshift(listselect(ivnotsimple),specorder-jj),planedata.planeeqs);
                end
                newIRcoordstoadd = zeros(length(listselect),3);
                newIRcoordstoadd(ivsimple,:) = newIRcoordssimple;
                newIRcoordstoadd(ivnotsimple,:) = newIRcoordsnotsimple;
                bigIRcoords(listselect,:) = newIRcoordstoadd;
	
            end
		
		end
        clear nspecular_select
    else
        nspecular_select_not_empty = 0;
    end
end

%   #######################################################################
%   #######################################################################
%   #######################################################################
%   #######################################################################
%   #######################################################################

%   ##################################################
%   #                                                #
%   #  S - (spec) - edge - spec - spec - R cases     #
%   #                                                #
%   #   (Pre- and) Postspec cases                    #
%   #                                                #
%   ##################################################
%
% Possible edges for S-E-spec-spec-R are seen (at least partly) by the receiver.
%
% The visibility must be checked, and also obstructions, but only
% for the paths between the edge and the receiver

% The vector masterivlist will always refer to the original data vector,
% i.e. PotentialIR, OriginsfromIR, reflorderIR etc
%
% First we pick out those indices where there was a single diffraction, 
% which wasn't last in the sequence 

if specorder > 1 && nspecular_select_not_empty == 1

    % The ii-value refers to the number of specular reflections after the
    % diffraction

	for ii = 1:specorder-1
        
        if showtext >= 3
            disp(['      Checking for ',int2str(ii),' spec refl after the edge diff'])    
        end
	
        % The selection of cases with ii specular reflection after one
        % diffraction have already been done, in the IR-building process.
	
        eval(['masterivlist = masterivlistorig(iv',int2str(ii),');'])
        
        if nedges < 256
            possibleedges = uint8(double(PotentialISESshift(masterivlist,specorder-ii)) - nplanes);
        elseif nedges < 65536
            possibleedges = uint16(double(PotentialISESshift(masterivlist,specorder-ii)) - nplanes);
        else
            possibleedges = uint32(double(PotentialISESshift(masterivlist,specorder-ii)) - nplanes);
        end
        
        % NB! possiblecombs will contain the specular sequences after the
        % diffraction but they appear in reversed order!
        possiblecombs = PotentialISESshift(masterivlist,specorder-ii+1:specorder);   

        % Specular combinations before the diffraction
        nmaxprespecs = specorder-1-ii;
        if nmaxprespecs == 0
            if nedges+nplanes < 256
                prespeccombs = uint8(zeros(length(masterivlist),1));
            elseif nedges+nplanes < 65536
                prespeccombs = uint16(zeros(length(masterivlist),1));
            else
                prespeccombs = uint32(zeros(length(masterivlist),1));                
            end
            nmaxprespecs = 1;
        else
            prespeccombs = PotentialISESshift(masterivlist,1:nmaxprespecs);
        end
       % Visibility of the edge in each sequence, seen from the source
        edgeweightsfromsou = ISESVISIBILITY(masterivlist);
        
        %--------------------------------------------------------------
        % Expand to take the edge subdivisions into account - 
        % treat all the edge subdivisions as separate edges for the
        % visibility and obstruction test.
        
        nposs = length(masterivlist);
        nposs = nposs*nedgesubs;
       
        expandedposscombs = repmat(possiblecombs,[1,nedgesubs]);
        expandedposscombs = reshape(expandedposscombs.',ii,nposs).';
  
        expandedprespecs = repmat(prespeccombs,[1,nedgesubs]);
        expandedprespecs = reshape(expandedprespecs.',nmaxprespecs,nposs).';
                        
        expandededgeweightsfromsou = repmat(edgeweightsfromsou,[1,nedgesubs]);
        expandededgeweightsfromsou = reshape(expandededgeweightsfromsou.',1,nposs).';
        
        expandedpossedges = repmat(possibleedges,[1,nedgesubs]);
        expandedpossedges = reshape(expandedpossedges.',1,nposs).';

        expandedmasterivlist = repmat(masterivlist,[1,nedgesubs]);
        expandedmasterivlist = reshape(expandedmasterivlist.',1,nposs).';

        if showtext >= 3
			disp(['         ',int2str(nposs),' Edge+IR segments found initially '])
        end
	        
        % Go through, iteratively, and check if the path from edge to R passes
        % through all reflection planes along the way
    
       for jj = ii:-1:1
	
            if nposs > 0
	
                if jj == ii
                    fromcoords = full(bigIRcoords(expandedmasterivlist,:));
                    [toedgecoords,edgeweightlist,edgenumberlist] = EDgetedgepoints(edgedata.edgestartcoords(possibleedges,:),edgedata.edgeendcoords(possibleedges,:),edgedata.edgelengthvec(possibleedges,:),nedgesubs);
                    tocoords = toedgecoords;
                else
                    eval(['tocoords = reflpoints',int2str(jj+1),';'])    
                    startlistvec = (1:length(expandedmasterivlist));                    
                    ivref = IRoriginsfrom(expandedmasterivlist);
                    for kk = jj:ii-2
                        ivnoIRexist = find(ivref==0);
                        if ~isempty(ivnoIRexist)
                            ivIRexist = startlistvec;
                            ivIRexist(ivnoIRexist) = [];
                            ivref(ivIRexist) = IRoriginsfrom(ivref(ivIRexist));
                        else
                            ivref = IRoriginsfrom(ivref);                            
                        end
                    end
                    ivnoIRexist = find(ivref==0);                   
                    if isempty(ivnoIRexist)
                        fromcoords = full(bigIRcoords(ivref,:)); 
                    else
                        ivIRexist = startlistvec;
                        ivIRexist(ivnoIRexist) = [];
                        fromcoords = zeros(length(ivref),3);
                        fromcoords(ivIRexist,:) = full(bigIRcoords(ivref(ivIRexist),:));

                        ISEScutout = PotentialISESshift(expandedmasterivlist(ivnoIRexist),specorder-ii+1:specorder);
                        newIRcoords = EDfindis(R,ISEScutout(:,ii),planedata.planeeqs);
                        for kk = 2:jj
                            newIRcoords = EDfindis(newIRcoords,ISEScutout(:,ii-kk+1),planedata.planeeqs);
                        end
                        fromcoords(ivnoIRexist,:) = newIRcoords;
                        
                    end
                    
                 end
                
                colno = ii-jj+1;
	
                [hitplanes,reflpoints,edgehits,edgehitpoints,edgehitnumbers,cornerhits,cornerhitpoints,cornerhitnumbers] = EDchkISvisible(fromcoords,tocoords,planedata.planeeqs(expandedposscombs(:,colno),4),planenvecs(expandedposscombs(:,colno),:),planedata.minvals(expandedposscombs(:,colno),:),...
				    planedata.maxvals(expandedposscombs(:,colno),:),planedata.planecorners(expandedposscombs(:,colno),:),planedata.corners,planedata.ncornersperplanevec(expandedposscombs(:,colno)));
                if ~isempty(edgehits) || ~isempty(cornerhits)
                    disp('WARNING! An edgehit or cornerhit occurred during the visibility test but this is not')
                    disp('         handled correctly yet.')
                end
                eval(['reflpoints',int2str(jj),' = reflpoints;'])
              
                expandedmasterivlist  = expandedmasterivlist(hitplanes);
                expandedposscombs     = expandedposscombs(hitplanes,:);
	%%            expandedreftoIRcoords = expandedreftoIRcoords(hitplanes);
                expandedpossedges     = expandedpossedges(hitplanes);
                expandedprespecs      = expandedprespecs(hitplanes,:);
                edgeweightlist        = edgeweightlist(hitplanes);
                expandededgeweightsfromsou = expandededgeweightsfromsou(hitplanes); 
                toedgecoords = toedgecoords(hitplanes,:);
%%%                tocoords          = tocoords(hitplanes,:);
	
                if jj < ii
                    for kk = jj+1:ii
                       eval(['reflpoints',int2str(kk),' = reflpoints',int2str(kk),'(hitplanes,:);'])
                    end
                end
	
                nposs = length(expandedmasterivlist);
                
            end
            if showtext >= 3
			    disp(['         ',int2str(nposs),' Edge+IR segments survived the visibility test in refl plane ',int2str(jj)])               
                if showtext >= 5
                   POTENTIALISES(expandedmasterivlist,:)
                end
            end
            
        end
        
        if obstructtestneeded && nposs > 0
        
            % Check obstructions for all the paths: edge -> plane1 -> plane2 -> ...
            % -> edge   (S to edge is not needed!)
            for jj = 1:ii+1
                if nposs > 0
                    if jj==1
                        fromcoords = R;    
                        startplanes = [];    
                    else
                        startplanes = expandedposscombs(:,ii-jj+2);
                        eval(['fromcoords = reflpoints',int2str(jj-1),';'])
                    end
                    if jj == ii+1
                        tocoords = toedgecoords;
                        endplanes = [edgedata.planesatedge(expandedpossedges,1) edgedata.planesatedge(expandedpossedges,2)];
%                        endplanes = [];    
                    else
                        eval(['tocoords = reflpoints',int2str(jj),';'])    
                        endplanes = expandedposscombs(:,ii-jj+1);    
                    end
                                        
                    [nonobstructedpaths,nobstructions,edgehits,cornerhits] = EDcheckobstrpaths(fromcoords,tocoords,startplanes,endplanes,planedata.canplaneobstruct,planedata.planeseesplane,...
                        planedata.planeeqs,planenvecs,planedata.minvals,planedata.maxvals,planedata.planecorners,planedata.corners,planedata.ncornersperplanevec,planedata.rearsideplane);
                    if ~isempty(edgehits) || ~isempty(cornerhits)
                        disp('WARNING! An edgehit or cornerhit occurred during the obstruction test but this is not')
                        disp('         handled correctly yet.')
                    end
	
                    if nobstructions > 0
                        expandedmasterivlist  = expandedmasterivlist(nonobstructedpaths);
                        expandedposscombs     = expandedposscombs(nonobstructedpaths,:);
                        expandedpossedges     = expandedpossedges(nonobstructedpaths);
                        expandedprespecs      = expandedprespecs(nonobstructedpaths,:);
                        edgeweightlist        = edgeweightlist(nonobstructedpaths);
                        expandededgeweightsfromsou = expandededgeweightsfromsou(nonobstructedpaths); 
%                        toedgecoords          = tocoords(nonobstructedpaths,:);
                        toedgecoords          = toedgecoords(nonobstructedpaths,:);
                        nposs = length(expandedmasterivlist);
                        
                        for kk = 1:ii
                            eval(['reflpoints',int2str(kk),' = reflpoints',int2str(kk),'(nonobstructedpaths,:);'])    
                        end
                        
                    end
	
                end
                if showtext >= 3
			       disp(['         ',int2str(nposs),' Edge+IR segments survived the obstruction test for path segment ',int2str(jj)])
                    if showtext >= 5
                        POTENTIALISES(expandedmasterivlist,:)
                    end
                end
           end      
 
        end
                
    	% % %     % Now we could in principle remove all the combinations we have been working on from the list masterivlistorig
	    % % %     eval(['masterivlistorig(iv',JJ(ii,:),',:) = [];'])
	           
        if nposs > 0

            % Visibility of edge segments is given by bitand-ing the
            % edge visibility from the (image) source and from the receiver

            edgeweightlist = bitand(edgeweightlist,expandededgeweightsfromsou);
            iv = find(edgeweightlist==0);

            if ~isempty(iv)
                expandedpossedges(iv) = [];
                expandedposscombs(iv,:) = [];
                expandedprespecs(iv,:) = [];
                edgeweightlist(iv) = [];
                expandedmasterivlist(iv) = [];
                nposs = length(edgeweightlist);
                
            end
             
            % Add the newly found combinations to the outdata list
            
            nnew = size(expandedprespecs,1);
            if nnew == 1
                colstoclear = find( expandedprespecs==0 );
            else
                colstoclear = find(sum(expandedprespecs>0)==0);
            end
            if ~isempty(colstoclear)
                expandedprespecs(:,colstoclear) = [];    
            end

            ncolsexist = size(prespeclist,2);
            ncolsnew = size(expandedprespecs,2);
            ncolstoadd = ncolsexist - ncolsnew;
            if ncolstoadd > 0
                prespeclist =  [prespeclist; [expandedprespecs zeros(nposs,ncolstoadd)]];
            else
                prespeclist =  [prespeclist; expandedprespecs];                
            end
            edgedifflist = [edgedifflist;expandedpossedges];
            postspeclist = [postspeclist;[expandedposscombs zeros(nposs,specorder-1-ii)]];
            bigedgeweightlist = [bigedgeweightlist;edgeweightlist];
            validIRcoords = [validIRcoords;bigIRcoords(expandedmasterivlist,:)];
            % NB! It is correct below that the indices for the ISCOORDS should be
            % ORIGINSFROM(masterivlist), rather than masterivlist.
            % The combinations in POTENTIALISES(masterivlist,:) all have
            % spec-spec-...-diff combinations and then
            % ISCOORDS(masterivlist,:) are zeros since a comb. that
            % ends with a diff has no image source. 
            % If we have a spec-spec-...-diff-spec comb., then we must
            % use ORIGINSFROM(ORIGINSFROM(masterivlist)) etc.
            ivref = ORIGINSFROM(expandedmasterivlist);
            for kk = 1:ii
                ivref = ORIGINSFROM(ivref);
            end
            validIScoords = [validIScoords;ISCOORDS(ivref,:)];

        end
	end

    clear PotentialISESshift bigIRcoords
end

%   #######################################################################
%   #
%   #   Pack the edge segments together because every little edge segment
%   #   is present as a separate edge
%   #   This can be done for all combinations at once.
%   #
%   #######################################################################

test = [prespeclist edgedifflist postspeclist];

if ~isempty(test)

	ncombs = length(edgedifflist);
	dtest = diff([0;prod(   double(test(1:ncombs-1,:)==test(2:ncombs,:)).'  ).']);
	ivremove = find(dtest==1);
	
	while ~isempty(ivremove)
        bigedgeweightlist(ivremove+1) = double(bigedgeweightlist(ivremove+1)) + double(bigedgeweightlist(ivremove));
        bigedgeweightlist(ivremove) = [];
        edgedifflist(ivremove) = [];
        prespeclist(ivremove,:) = [];
        postspeclist(ivremove,:) = [];
        validIScoords(ivremove,:) = [];
        validIRcoords(ivremove,:) = [];
	
        test = [prespeclist edgedifflist postspeclist];
        ncombs = length(edgedifflist);
        dtest = diff([0;prod(   double(test(1:ncombs-1,:)==test(2:ncombs,:)).'  ).']);
        ivremove = find(dtest==1);
	end
end

%   #######################################################################
%   #
%   #   The weights of the visible edge segments should be
%   #   translated into start and end points, together with the visibility
%   #   weights from the receiver.
%   #   This can be done for all combinations at once.
%   #
%   #######################################################################

% As a start we set the start and end values to 0 and 1, i.e. assuming full
% visibility.

if specorder >= 1 && ~isempty(edgedifflist)
	startandendpoints = [startandendpoints;[zeros(length(edgedifflist),1) ones(length(edgedifflist),1)]];
	
    totvisibility = bigedgeweightlist;

	ivnovisibility = find(totvisibility==0);
	if ~isempty(ivnovisibility)
        disp('Clearing some cases.')
        edgedifflist(ivnovisibility) = [];
        prespeclist(ivnovisibility,:) = [];
        postspeclist(ivnovisibility,:) = [];
        bigedgeweightlist(ivnovisibility) = [];
        totvisibility(ivnovisibility) = [];
        validIScoords(ivnovisibility,:) = [];
        validIRcoords(ivnovisibility,:) = [];
	end    
	
	ivtoaccountfor = [1:length(totvisibility)].';
	
	% If there are edges with full visibility we don't need to do anything to
	% the start and end values, but we remove the edges from the list to
	% process.
	
	ivwholeedges = find( totvisibility == maxvisibilityvalue);
	if ~isempty(ivwholeedges)
        ivtoaccountfor(ivwholeedges) = [];
	end
	
	if ~isempty(ivtoaccountfor)
        ncombs = length(ivtoaccountfor);
        bitpattern = zeros(ncombs,nedgesubs);
        for ii=1:nedgesubs
            bitpattern(:,ii) = bitget(totvisibility(ivtoaccountfor),ii); 
        end
        dbit1 = diff([zeros(ncombs,1) bitpattern].').';
        dbit2 = [dbit1 -bitpattern(:,nedgesubs)]; 
        
        nsegments = ceil((sum(abs(dbit1.')).')/2);
	
        ivonesegments = find(nsegments==1);
        if ~isempty(ivonesegments)

            nonesegments = length(ivonesegments);
            multvec = 2.^(0:nedgesubs);
            segstartpos = round(log(sum( ((dbit2(ivonesegments,:)== 1).*multvec(ones(nonesegments,1),:)).').')/log(2))+1;
            segendpos   = round(log(sum( ((dbit2(ivonesegments,:)==-1).*multvec(ones(nonesegments,1),:)).').')/log(2))+1;
	
            ivmodify = find(segstartpos==1);
            segstartpos(ivmodify) = ones(size(ivmodify))*1.5;
            ivmodify = find(segendpos>nedgesubs);
            segendpos(ivmodify) = ones(size(ivmodify))*(nedgesubs+0.5);
	
            startandendpoints(ivtoaccountfor(ivonesegments),1) = (segstartpos-1.5)/(nedgesubs-1);
            startandendpoints(ivtoaccountfor(ivonesegments),2) = (segendpos-1.5)/(nedgesubs-1);
                    
        end    
        
        % If we have some two-or-more-subsegments cases, they will be
        % discovered by the if-condition below
        
        if length(ivonesegments) < ncombs
            for nsegmentstocheck = 2:ceil(nedgesubs/2)
%                disp(['Checking for ',int2str(nsegmentstocheck),' sub-segments'])
     
                ivNsegments = find(nsegments==nsegmentstocheck);
                if ~isempty(ivNsegments)
                    [n1,n2] = size(startandendpoints);
                    if n2 < 2*nsegmentstocheck
                        startandendpoints = [startandendpoints zeros(n1,2*nsegmentstocheck-n2)];    
                    end
                    for jj = 1:length(ivNsegments)
                        ivstartbits = find(dbit2(ivNsegments(jj),:) == 1);
                        ivstartbits = (ivstartbits==1)*1.5 + (ivstartbits~=1).*ivstartbits;
                        ivendbits = find(dbit2(ivNsegments(jj),:) == -1);
                        ivendbits = (ivendbits>nedgesubs)*(nedgesubs+0.5) + (ivendbits<=nedgesubs).*ivendbits;
                        
                        for kk = 1:nsegmentstocheck                                                        
                            startandendpoints(ivtoaccountfor(ivNsegments(jj)),(kk-1)*2+1) = (ivstartbits(kk)-1.5)/(nedgesubs-1);
                            startandendpoints(ivtoaccountfor(ivNsegments(jj)),(kk-1)*2+2) = (ivendbits(kk)-1.5)/(nedgesubs-1);
                        end
                    end                
                end
                
            end
        end        
	end
end

%   #######################################################################
%   #
%   #   Construct a list guide, which will tell which rows have only
%   #   d, which rows have sd etc
%   #   
%   #######################################################################

[n1,n2] = size(prespeclist);
if n1 ~= 0
	if n2 == 0
        prespeclist = zeros(n1,1);    
	end
	
	if n2 > 1
        nprespecs  = sum(prespeclist.' > 0).';
	else
        nprespecs = (prespeclist>0);
	end
	n2 = size(postspeclist,2);
	if n2 > 1
        npostspecs  = sum(postspeclist.' > 0).';
	else
        npostspecs = (postspeclist>0);    
	end
	
	listofdiffcol   = 1 + nprespecs;
	listofreflorder = listofdiffcol + npostspecs;
	
	[B,ivec,jvec] = unique([listofreflorder listofdiffcol],'rows');
	nuniquecombs = length(ivec);
	ntotcombs = length(jvec);
	
	listguide = zeros(nuniquecombs,3);
	listoforders = zeros(nuniquecombs,2);
	sortvec = zeros(ntotcombs,1);
	for ii = 1:length(ivec)
        ivfindcombs = find(jvec==ii);
        listguide(ii,1) = length(ivfindcombs);
        if ii > 1
            listguide(ii,2) = listguide(ii-1,3)+1;
        else
            listguide(ii,2) = 1;    
        end
        listguide(ii,3) = listguide(ii,2)+listguide(ii,1)-1;
        listoforders(ii,:) = [B(ii,1) B(ii,2)]; 
	
        sortvec(listguide(ii,2):listguide(ii,3)) = ivfindcombs;
        
	end
	
	prespeclist = prespeclist(sortvec,:);
	postspeclist = postspeclist(sortvec,:);
	validIScoords = validIScoords(sortvec,:);
	validIRcoords = validIRcoords(sortvec,:);
	edgedifflist = edgedifflist(sortvec,:);
	startandendpoints = startandendpoints(sortvec,:);
	bigedgeweightlist = bigedgeweightlist(sortvec,:);
    
    % The prespeclist needs to be shifted back to the left because it is
    % "right-aligned"
    
    [ncombs,ncols] = size(prespeclist);
    if ncols > 1
        if ncombs > 1
            iv = find( ( prespeclist(:,1)==0 ) & ( prespeclist(:,2)~=0 ) );
            if ~isempty(iv)
               prespeclist(iv,1:ncols-1) = prespeclist(iv,2:ncols);
               prespeclist(iv,ncols) = 0;
            end

            for jj = 3:ncols
                iv = find( ( sum(prespeclist(:,1:jj-1).').' == 0 ) & ( prespeclist(:,jj)~=0 ) );
                if ~isempty(iv)
                   prespeclist(iv,1:ncols-(jj-1)) = prespeclist(iv,jj:ncols);
                   prespeclist(iv,ncols-(jj-2):ncols) = 0;
                end
            end    
        else
           if prespeclist(1,1) == 0               
                ninitialzeros =  find(diff(cumsum(prespeclist)));
                if ~isempty(ninitialzeros)
                    ninitialzeros = ninitialzeros(1);
                    prespeclist(1:ncols-ninitialzeros) = prespeclist(ninitialzeros+1:ncols);
                    prespeclist(ncols-ninitialzeros+1:ncols) = 0;
                end
           end
        end
    end
    
else
    prespeclist = [];
    postspeclist = [];
    validIScoords = [];
    validIRcoords = [];
    edgedifflist = [];
    startandendpoints = [];
    bigedgeweightlist = [];
    listguide = [];
    
end
