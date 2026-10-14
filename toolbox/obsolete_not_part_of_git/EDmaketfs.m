function [tfdirect,tfgeom,tfdiff] = EDmaketfs(envdata,planedata,edgedata,edgetoedgedata,...
    pathstruct,specorder,frequencies,Rstart,guiderowstouse,includedirectsound,showtext)
% EDmaketfs - Constructs transfer functions from a list of paths.
%
% Input parameters:
%   envdata,planedata, edgedata, edgetoedgedata           structs
%   pathstruct          struct with these fields:
%       pathtypevec
%       reflpaths
%       specextradata
%       edgeextradata
%       S
%       R
%       mainlistguide
%       mainlistguidepattern
%       directsoundrow
%       allspecrows
%       firstdiffrow
%       Sinsideplanenumber
%       Rinsideplanenumber
% 	specorder           The maximum order of reflections that is calculated.
%   frequencies         A list of frequencies that the TF should be
%                       computed for.
% 	Rstart              The reference distance that the distance zero of
%                       the transfer response refers to.
%   guiderowstouse      (optional) A list of values, indicating which rows in the
%                       mainlistguide and mainlistguidepattern that should
%                       be used. This way it is possible to select only
%                       diffraction, for instance. If this list is not
%                       specified, all components will be used.
%   includedirectsound  (optional) 0 or 1, indicating whether the direct
%                       sound should be included or not. Default: 1 
%   CAIR, showtext   Global parameters.
%
% Output parameters:
%   tfdirect            The TF containing the direct sound, size [nfrequencies,1].
%   tfgeom              The TF containing the specular reflections, size [nfrequencies,1].
%   tfdiff              The TF containing the diffracted components, size [nfrequencies,1].
%
% Uses functions EDcalcdist EDwedge1st_fd 
% 
% Peter Svensson 29 Nov. 2017 (peter.svensson@ntnu.no)
%
% [tfdirect,tfgeom,tfdiff] = EDmaketfs(envdata,planedata,edgedata,edgetoedgedata,...
%     reflpaths,specextradata,edgeextradata,S,R,mainlistguide,mainlistguidepattern,...
%     directsoundrow,allspecrows,firstdiffrow,Sinsideplanenumber,Rinsideplanenumber,specorder,...
%     frequencies,Rstart,guiderowstouse,includedirectsound,showtext);

% 16 Jul 2010   Functioning version
% 3 Dec 2014    Changed to using structs for input
% 31 March 2015 Cleaned up the help text and removed some redundant input
%               parameters.
% 1 Nov. 2016   Had forgotten one change of planesatedge to edgedata.planesatedge
% 2 March 2017 Checked against ESIE0maketfs
% 21 Oct. 2017 Bug with boost factor!
% 27 Nov. 2017 Copied from ESIE2toolbox
% 28 Nov. 2017 Cleaned up code. Introduced the non-glbal input parameter
%               showtext. Introduced the input struct "pathstruct"
% 29 Nov. 2017 Changed two calls to ESIE2coordtrans

reflpaths = pathstruct.reflpaths;
specextradata = pathstruct.specextradata;
edgeextradata = pathstruct.edgeextradata;
S = pathstruct.S;
R = pathstruct.R;
mainlistguide = pathstruct.mainlistguide;
mainlistguidepattern = pathstruct.mainlistguidepattern;
directsoundrow = pathstruct.directsoundrow;
allspecrows = pathstruct.allspecrows;
firstdiffrow = pathstruct.firstdiffrow;
Sinsideplanenumber = pathstruct.Sinsideplanenumber; 
Rinsideplanenumber = pathstruct.Rinsideplanenumber;

%--------------------------------------------------------------------------
% We look for the optional vector multfactors .
% If it was there, then those values will be used to switch off, or boost,
% any reflection path.
% Those values will only be used for diffraction paths.

if exist('multfactors','var') == 0
    multfactors = ones(size(reflpaths,1),1);    
else
    if size(multfactors,1) ~= size(reflpaths,1)
        error(['The edpaths file contained a vector multfactors which does not have the same length as reflpaths.'])    
    end
end

%--------------------------------------------------------------------------

% Sdirection = [1 0 0];
% maxnedgesorplanes = max(max(reflpaths));
% if ~isempty(maxnedgesorplanes)
% %     multfac = 10^ceil(log10(double(maxnedgesorplanes)+1));
% else
% %     multfac = 0;
% end

if nargin < 11
    showtext = 0;
    if nargin < 10
        includedirectsound = 1;
        if nargin < 9
            guiderowstouse = [];
        end
    end
end

if isempty(guiderowstouse)
    usesubset = 0;    
    if showtext > 2
        disp('      Using all components')
    end
else
    usesubset = 1;        
    if showtext > 2
        disp('      Using only some components')
    end
end

nyvec = pi./(2*pi - edgedata.closwedangvec);    
onesvec1 = ones(1,3);

souspecboost = 1;
if ~isempty(Sinsideplanenumber)
    if reflfactors(Sinsideplanenumber(1)) ~= 0
        souspecboost = 2;
    end
end

recspecboost = 1;
if ~isempty(Rinsideplanenumber)
    if reflfactors(Rinsideplanenumber(1)) ~= 0
        recspecboost = 2;
    end
end

if isstruct(edgetoedgedata) == 0
    userwantsdiff2 = 0;
else
    if isempty(edgetoedgedata.re1sho)
        userwantsdiff2 = 0;
    else
        userwantsdiff2 = 1;    
    end
end

nfrequencies = length(frequencies);

%-------------------------------------------------------

if exist('mainlistguide','var') ~= 1
    mainlistguide = [];    
else
    mainlistguide = double(mainlistguide);    
end

[ncomponents,ncols] = size(reflpaths);
nrowsinlist = size(mainlistguide,1);
if ncols > 1
    difforderinlist = sum(mainlistguidepattern.'=='d').';
else
    difforderinlist = (mainlistguidepattern=='d');    
end

lastNdiffrow = zeros(1,specorder);
for ii = 1:specorder
    iv = find(difforderinlist <= ii );
    if ~isempty(iv)
        lastNdiffrow(ii) = iv(end);
    end
end

% nrefltogetherwdiff = sum(mainlistguidepattern.'=='d').'.*(  sum(mainlistguidepattern.'=='s').');

if showtext >= 3
	disp(['   Constructing TF. ',int2str(ncomponents),' components.'])
end

tfdirect = zeros(nfrequencies,1);
tfgeom   = zeros(nfrequencies,1);
tfdiff   = zeros(nfrequencies,1);
% tftot    = zeros(nfrequencies,1);

if ncomponents == 0
    return
end

%##############################################################################
%##############################################################################
%##############################################################################
%##############################################################################
%
%  Diffraction once, possibly with pre- and post-specular reflections
%
%##############################################################################

directsoundonboundary = 0;
specreflonboundary = zeros(size(planedata.planeeqs,1),1);

if firstdiffrow ~= 0
    
    % First we remove the 'd' for all the rows in the mainlistguidepattern that
    % we do not want to use.
    
    if usesubset == 1
        singdiffrows = find(difforderinlist==1);
        rowstoclear = setdiff(singdiffrows,guiderowstouse);
        if ~isempty(rowstoclear)
            mainlistguidepattern(rowstoclear,:) = mainlistguidepattern(rowstoclear,:)*0;    
        end
    end

    % Here we can use the single-diffraction calculation for all
    % combinations. The source should either be the original source or
    % an image source. The receiver should either be the original
    % receiver or an image receiver.
    
    % First we create lists that specify whether there are specular
    % reflections before the edge diffraction, and after the edge
    % diffraction. Also, a list which gives the column where the edge
    % diffraction is present can be used to efficiently extract the
    % edge number. 
    % NB! The lists prespeclist, postspeclist and singdiffcol all have
    % the (short) length of the mainlistguide. 
    % The data in these lists will then be used to create the long
    % lists prespecpresent and postspecpresent which simply are 0 or 1.
    % A matrix called edgemask will have the same number of columns as
    % the reflpaths matrix. For each row, there will be a 1 in the
    % column where the edge diffraction is so we can efficiently find
    % the edge number.

    multmatrix = (1:specorder);
    if size(mainlistguidepattern,2) > specorder
        multmatrix = (1:size(mainlistguidepattern,2));
    end
    multmatrix = multmatrix( ones(nrowsinlist,1),:  );
    
    if ncols > 1
        singdiffcol = sum( (multmatrix.*(mainlistguidepattern=='d')).' ).';
        singdiffcol = singdiffcol.*(sum(mainlistguidepattern.'=='d').'<=1);
        nrefl = sum( (mainlistguidepattern=='s' | mainlistguidepattern=='d').' ).';
    else
        singdiffcol = mainlistguidepattern=='d';
        nrefl = singdiffcol | mainlistguidepattern=='s';
    end
    
    prespeclist = singdiffcol-1;
    prespeclist = prespeclist.*(prespeclist>0);
    postspeclist = (nrefl-singdiffcol).*(singdiffcol~=0);
    
    prespecpresent = zeros(ncomponents,1);
    postspecpresent = zeros(ncomponents,1);
    edgemask = zeros(ncomponents,specorder);

    diffrowsthatareOK = (firstdiffrow:lastNdiffrow(1));
    if usesubset == 1
        diffrowsthatareOK = intersect(diffrowsthatareOK,guiderowstouse);
    end
    for ii = diffrowsthatareOK    %firstdiffrow:lastNdiffrow(1),
        iv = ((mainlistguide(ii,2)):(mainlistguide(ii,3)));
        onesvec2 = ones(length(iv),1);
        edgemask(iv,singdiffcol(ii)) = onesvec2;            
%         if ii > specorder+2,
            prespecpresent(iv) = (prespeclist(ii)>0)*onesvec2;
            postspecpresent(iv) = (postspeclist(ii)>0)*onesvec2;
%         end
    end
    prespecpresent = prespecpresent(:,onesvec1);
    postspecpresent = postspecpresent(:,onesvec1);

    if ncols > 1
%         edgenumb = sum( (edgemask.*double(reflpaths)).').';
        edgenumb = sum( (edgemask.*double(reflpaths(:,1:specorder))).').';
    else
        edgenumb = edgemask.*double(reflpaths);    
    end

%     nedgetfs = size(edgenumb,1);
%     nnonzeroedgetfs = sum(edgenumb(:,1)>0);

%     tfcounter = 0;
    edgelist = unique(edgenumb);
    for ii = 1: length(edgelist)
        edge = edgelist(ii);
        if edge~= 0
           iv = find(edgenumb==edge);
           ncombs = length(iv);
           
           if ncombs > 0 && any(multfactors(iv))
               if showtext >= 4
                    disp(['   Edge ',int2str(edge),': ',int2str(ncombs),' combinations'])    
               end               
                 
%                 onesvec3 = ones(ncombs,1);
	
                IS = full(specextradata(iv,1:3));
                IR = full(specextradata(iv,4:6));                
            
                edgestart = full(edgeextradata(iv,1));
                edgeend   = full(edgeextradata(iv,2));
                
                % Calculate the cyl coordinates of all IS
                % 
                % Calculate the cyl coord of all IR
	
                edgecoords = [edgedata.edgestartcoords(edge,:);edgedata.edgeendcoords(edge,:)];
            	[rs,thetas,zs,rr,thetar,zr] = EDcoordtrans2(IS,IR,edgecoords,edgedata.edgenvecs(edge,:));
	
				bc = planedata.reflfactors(edgedata.planesatedge(edge,:)).';

                for jj = 1:ncombs
                    if multfactors(iv(jj)) > 0

                       [tfnew,singularterm] = EDwedge1st_fd(envdata.cair,frequencies,edgedata.closwedangvec(edge),rs(jj),thetas(jj),zs(jj),rr(jj),thetar(jj),zr(jj),...
                         edgedata.edgelengthvec(edge)*[edgestart(jj) edgeend(jj)],'n',Rstart,bc);                  

                          if any(singularterm)
                              
                              if singularterm(2) || singularterm(3)
                                  directsoundonboundary = 1;
                              elseif singularterm(1)
                                if specorder == 1
                                   specreflonboundary(edgedata.planesatedge(edge,2)) = 1;
                               else
                                   disp('WARNING! A specular refl. of order > 1 is half obscured but this is not handled yet!');    
                                end
                              elseif singularterm(4)
                                  if specorder == 1
                                       specreflonboundary(edgedata.planesatedge(edge,1)) = 1;
                                   else
                                       disp('WARNING! A specular refl. of order > 1 is half obscured but this is not handled yet!');    
                                  end
                              end
                          end
		
                    
                        % Decide whether the TF should be boosted or not
                        %
                        % The factors souspecboost and recspecboost have
                        % values other than 1 if the source or receiver is directly
                        % at a plane, for the boosting of the direct sound, and
                        % specular reflections.
                        %
                        % This boost factor should be used also for edge
                        % diffraction if:
                        %   1.  There are some specular reflections between the 
                        %       source and the edge
                        %   or
                        %   2.  There are no specular reflections between the
                        %       source and the edge, but the source/receiver is at
                        %       a plane which does not belong to the edge.
                        
                        if souspecboost ~= 1 || recspecboost ~= 1
                        
                            boostfactor = 1;
                            if prespecpresent(iv(jj)) == 1
                                boostfactor = souspecboost;    
                            else
                                if ~isempty(Sinsideplanenumber)
                                    if Sinsideplanenumber(1)~=planesatedge(edge,1) && Sinsideplanenumber(1)~=planesatedge(edge,2) 
                                        boostfactor = souspecboost;
                                    end
                                end
                            end
                            if postspecpresent(iv(jj)) == 1
                                boostfactor = boostfactor*recspecboost;    
                            else
                                if ~isempty(Rinsideplanenumber)
                                    if Rinsideplanenumber(1)~=planesatedge(edge,1) && Rinsideplanenumber(1)~=planesatedge(edge,2) 
                                        boostfactor = boostfactor*recspecboost;
                                    end
                                end
                            end
                            if boostfactor ~= 1
                                tfnew = tfnew*boostfactor;    
                            end
                        
                        end    % ....    if souspecboost ~= 1 | recspecboost ~= 1,


                        tfdiff = tfdiff + tfnew*double(multfactors(iv(jj)));



                    end   %  ...  if multfactors(iv) > 0,
                    
                end  % .....                for jj = 1:ncombs,

            end   %....    if ncombs > 0 & any(multfactors(iv)),
            
        end    % ....         if edge~= 0,
        
    end  % ....     for ii = 1: length(edgelist),
    
end   % ...if firstdiffrow ~= 0,

nspecreflonboundary = sum(specreflonboundary);


%##############################################################################
%##############################################################################
%##############################################################################
%##############################################################################
%
%      The direct sound
%
%##############################################################################

if usesubset == 1 && directsoundrow == 1
    if any(guiderowstouse == 1) == 0
        directsoundrow = 0;    
    end
end

if directsoundrow == 1 && includedirectsound == 1    
    dist = norm(R-S);
%     slotnumberfrac = (dist-Rstart)/CAIR*FSAMP+1;
    amp = souspecboost*recspecboost./dist;
    if directsoundonboundary
        amp = amp/2;    
        if showtext >= 4
            disp('HALVING DIRECT SOUND')
        end
    end
    tfdirect = amp*exp(-1i*2*pi*frequencies(:)/envdata.cair*(dist-Rstart));
    
	if ncomponents == 1
%         tftot = tftot + tfdirect;
        
        return
	end    
else
    if directsoundonboundary == 1
        
        dist = norm(R-S);
%         slotnumberfrac = (dist-Rstart)/CAIR*FSAMP+1;
        amp = souspecboost*recspecboost./dist;
        if directsoundonboundary
            amp = amp/2;    
            if showtext >= 4
                disp('HALVING DIRECT SOUND')
            end
        end
%         tfdirect = amp*exp(-j*2*pi*frequencies(:)/CAIR*(dist-Rstart));
        tfdirect = amp*exp(-1i*2*pi*frequencies(:)/envdata.cair*(dist-Rstart));
        
      
    end
    if showtext >= 4
        disp('      No direct sound')    
    end
end

%##############################################################################
%##############################################################################
%##############################################################################
%##############################################################################
%
%  All-specular			s, ss, sss, etc
%
%##############################################################################

if allspecrows(1) ~= 0
    if usesubset == 1
        specrowsthatareOK = intersect(allspecrows,guiderowstouse);
    else
        specrowsthatareOK = (allspecrows(1):allspecrows(2));
    end
    if ~isempty(specrowsthatareOK)
        if usesubset == 1
            temp = mainlistguide(specrowsthatareOK,2:3);
            ivspec = (double(mainlistguide(specrowsthatareOK(1),2)):double(mainlistguide(specrowsthatareOK(1),3)));
            for ii = 2:size(temp,1)
                ivspec = [ivspec  (double(mainlistguide(specrowsthatareOK(ii),2)):double(mainlistguide(specrowsthatareOK(ii),3)))];
            end
        else
            ivspec = (mainlistguide(allspecrows(1),2):mainlistguide(allspecrows(2),3));        
        end
        
        if nspecreflonboundary > 0
            listofspecreflonboundary = find(specreflonboundary);
            specrefltoadd = setdiff(listofspecreflonboundary,ivspec);
            if ~isempty(specrefltoadd)
                error('ERROR: We need to add some code here, to reinsert pruned spec. refl.')    
            end
        end
        
        if showtext >= 4
            disp(['      ',int2str(length(ivspec)),' specular reflections'])    
        end
%         dist = ESIE2calcdist(full(specextradata(ivspec,1:3)),R);
        dist = EDcalcdist(full(specextradata(ivspec,1:3)),R);
        if size(reflpaths,2) == 1
            refl = planedata.reflfactors(  reflpaths(ivspec) );
        else
        	if showtext >= 2
	           disp('WARNING: havent implemented reflfactors yet for spec reflorder >1!'); 
    	       disp('   reflfactors will be set to 1');
    	   end
           refl = 1;
        end
        
        nspec = length(dist);

        if nspecreflonboundary > 0
            specamp = ones(size(souspecboost));
            amplitudeshouldbehalf = ismember(reflpaths(ivspec,1),listofspecreflonboundary);
            specamp = specamp - amplitudeshouldbehalf*0.5;
            if nspec == 1
%                 tfnew = specamp*souspecboost*recspecboost*refl.*exp( -j*2*pi*frequencies(:)/CAIR*(dist-Rstart))/dist;		               
                tfnew = specamp*souspecboost*recspecboost*refl.*exp( -1i*2*pi*frequencies(:)/envdata.cair*(dist-Rstart))/dist;		               
            else
                dist = dist(:).';
                frequencies = frequencies(:);
%                tfnew = specamp*souspecboost*recspecboost*refl.*exp( -j*2*pi*frequencies(:,ones(1,nspec))/CAIR.*(dist(ones(nfrequencies,1),:)-Rstart))./dist(ones(nfrequencies,1),:);	 
               tfnew = specamp*souspecboost*recspecboost*refl.*exp( -1i*2*pi*frequencies(:,ones(1,nspec))/envdata.cair.*(dist(ones(nfrequencies,1),:)-Rstart))./dist(ones(nfrequencies,1),:);	 
               tfnew = sum(tfnew.').';
            end
        else
%     			irnew = souspecboost*recspecboost*ESIE2irfromslotlist((dist-Rstart)/CAIR*FSAMP+1,1./dist);		

            if nspec == 1
%                 tfnew = souspecboost*recspecboost*refl.*exp( -j*2*pi*frequencies(:)/CAIR*(dist-Rstart))/dist;		               
                tfnew = souspecboost*recspecboost*refl.*exp( -1i*2*pi*frequencies(:)/envdata.cair*(dist-Rstart))/dist;		               
            else
                dist = dist(:).';
                refl = refl(:).';
                frequencies = frequencies(:);
%                    tfnew = souspecboost*recspecboost*refl.*exp( -j*2*pi*frequencies(:,ones(1,nspec))/CAIR.*(dist(ones(nfrequencies,1),:)-Rstart))./dist(ones(nfrequencies,1),:);	 
% 				tfnew = souspecboost*recspecboost*refl(ones(nfrequencies,1),:).*exp( -j*2*pi*frequencies(:,ones(1,nspec))/CAIR.*(dist(ones(nfrequencies,1),:)-Rstart))./dist(ones(nfrequencies,1),:);	 
				tfnew = souspecboost*recspecboost*refl(ones(nfrequencies,1),:).*exp( -1i*2*pi*frequencies(:,ones(1,nspec))/envdata.cair.*(dist(ones(nfrequencies,1),:)-Rstart))./dist(ones(nfrequencies,1),:);	 
				tfnew = sum(tfnew.').';
            end
        end
        tfgeom = tfgeom + tfnew;
            
    else
        if nspecreflonboundary > 0
%             listofspecreflonboundary = find(specreflonboundary);
            error('ERROR: We need to add some code here, to reinsert pruned spec. refl., pos. 2')    
        end        
    end    
else
% % %     if nspecreflonboundary > 0,
% % %         listofspecreflonboundary = find(specreflonboundary);
% % %         
% % %         [xis] = ESIE2findis(S,listofspecreflonboundary,planeeqs,1,[1 1 1]);
% % %         
% % %         disp(['WARNING: We have some new code here, to reinsert pruned spec. refl., pos. 3'])
% % %         dist = ESIE2calcdist(xis,R);
% % %         irnew = 0.5*souspecboost*recspecboost*ESIE2irfromslotlist((dist-Rstart)/CAIR*FSAMP+1,1./dist);		                
% % %         ngeom = length(irgeom);
% % %         nnew = length(irnew);
% % %         if nnew > ngeom,
% % %            irgeom = [irgeom;zeros(nnew-ngeom,1)];
% % %         end
% % %         irgeom(1:nnew) = irgeom(1:nnew) + irnew;
% % %         
% % %     end        
end

%##############################################################################


%##############################################################################
%##############################################################################
%##############################################################################
%##############################################################################
%
%  Multiple diffraction, possibly with pre- and post-specular reflections
%
%##############################################################################

if userwantsdiff2 == 1
    
% % %     JJ = setstr(32*ones(specorder,1));
% % %     for jj=1:specorder
% % %         jjstr = int2str(jj);
% % %         JJ(jj,1:length(jjstr)) = jjstr;
% % %     end
% % %     [n1,n2] = size(JJ);

    edgeseespartialedge = edgetoedgedata.edgeseespartialedge;
    edgestartcoords = edgedata.edgestartcoords;
    edgeendcoords = edgedata.edgeendcoords;
    edgenvecs = edgedata.edgenvecs;
    reftoshortlistE = edgetoedgedata.reftoshortlistE;
    re1sho =edgetoedgedata.re1sho;
    re2sho = edgetoedgedata.re2sho;
    thetae1sho = edgetoedgedata.thetae1sho;
    thetae2sho = edgetoedgedata.thetae2sho;
    ze1sho = edgetoedgedata.ze1sho;
    ze2sho = edgetoedgedata.ze2sho;
    edgelengthvec = edgedata.edgelengthvec;
    
	for Ndifforder = 2:specorder
        if showtext >= 3
            disp(['   Diffraction order ',int2str(Ndifforder)])    
        end
	
        if any(difforderinlist==Ndifforder) % & elemsize(Ndifforder) > 0,
	
            % Calculate some general parameters that are shared for all
            % N-diffraction calculations
            
%             divmin = CAIR/(FSAMP*elemsize(Ndifforder));
%             ndivvec = ceil(abs( edgelengthvec.' )/divmin);
%             dzvec = (edgelengthvec.')./ndivvec;
	
            ncylrows = 4*(Ndifforder-1);
            
            % Here we can use the double-diffraction calculation for all
            % combinations. The source should either be the original source or
            % an image source. The receiver should either be the original
            % receiver or an image receiver.
            
            noffset = lastNdiffrow(Ndifforder-1);
            ivNdiff = (noffset+1:lastNdiffrow(Ndifforder));
            ndiffcombs = length(ivNdiff);
            zerosvec1 = zeros(ndiffcombs,1);
            
%             ndoubrows = length(ivNdiff);        
            previousrow = lastNdiffrow(Ndifforder-1);
            if previousrow > 0
                noffsetlonglist = mainlistguide(previousrow,3);
            else
                noffsetlonglist = 0;    
            end
            nremaining = ncomponents - (mainlistguide(lastNdiffrow(Ndifforder),3));
            nlonglist = ncomponents-noffsetlonglist-nremaining;
            ivlonglist = [mainlistguide(ivNdiff(1),2):mainlistguide(ivNdiff(end),3)].';
	
            diffcols = zerosvec1(:,ones(1,Ndifforder));
            for ii = ivNdiff(1):ivNdiff(end)
                diffcols(ii-ivNdiff(1)+1,:) = find(mainlistguidepattern(ii,:)=='d');
            end
            
            nreflorder = sum( (mainlistguidepattern(ivNdiff(1):ivNdiff(end),:)=='d'|mainlistguidepattern(ivNdiff(1):ivNdiff(end),:)=='s').').';
%             nprespecs = diffcols(:,1)-1;
            nmidspecs = diffcols(:,2)-diffcols(:,1)-1;
%             npostspecs = nreflorder-diffcols(:,Ndifforder);
	
            if ndiffcombs > 1
                ivreftolonglist = ivlonglist(:,ones(1,ndiffcombs));
                comppattern = mainlistguide(ivNdiff,2).';        
                comppattern = comppattern(ones(nlonglist,1),:);
            
                rowinpatternlist = sum( (ivreftolonglist>=comppattern).' ).';
            else
                rowinpatternlist = ones(size(ivlonglist));    
            end
            
            % Construct a long matrix with the edge numbers extracted from the
            % reflpaths matrix.
	
            longnmidspec = zeros(nlonglist,1);
            iv = find(nmidspecs>0);
            for ii = 1:length(iv)
                ivreftolonglist = (mainlistguide(ivNdiff(iv(ii)),2):mainlistguide(ivNdiff(iv(ii)),3)) - noffsetlonglist;
                longnmidspec(ivreftolonglist) = nmidspecs(iv(ii));
            end
            
            %------------------------------------------------------------------
            % edgepattern will be a matrix which, for each reflection combination, contains
            % the 2,3,4,... edge numbers that are involved in each path
            % regardless of there are specular reflections before, in the
            % middle, or after.
            
            edgepattern = zeros(nlonglist,Ndifforder);
            countvec = [1:nlonglist].';
            reflpathscut = reflpaths(ivlonglist,:);            
            for ii = 1:Ndifforder
                ivrefcol = countvec + (diffcols(rowinpatternlist,ii)-1)*nlonglist;
                edgepattern(:,ii) = reflpathscut(ivrefcol);
            end
	
            %------------------------------------------------------------------
            % midspecpattern will be a matrix which, for each reflection combination, contains
            % the 1,2,3,4,... specular reflection numbers (i.e., plane numbers) that are involved
            % in each path between diffraction 1 and diffraction 2.
            % First, the reflection component immediately after diffraction 1
            % is selected for every reflection. Then there is a masking which
            % zeroes the midspecpattern value for all the combinations that
            % don't have any midspec.
            % NB! For combinations like ssssdd, the first step will point to a
            % column which is outside the matrix. This is fixed by maximizing
            % the column number.
            
            midspecpattern = zeros(nlonglist,max([specorder-Ndifforder 1]));
            maxpointvalue = nlonglist*max([specorder-Ndifforder 1]);
            for ii = 1:specorder-Ndifforder
                ivrefcol = countvec + (diffcols(rowinpatternlist,1)+(ii-1))*nlonglist;
                ivrefcol = mod(ivrefcol-1,maxpointvalue)+1;
                midspecpattern(:,ii) = double(reflpathscut(ivrefcol)).*(longnmidspec>=ii);
            end
%             diff1col = diffcols(rowinpatternlist,1);
	
            
            % Many combinations include the same edge-pair/N-let, so we extract the
            % unique edge pairs/N-lets, and go through them in the ii-loop
            
            [edgeshortlist,ivec,jvec] = unique(edgepattern,'rows');
%             lastndivcomb = zeros(1,Ndifforder);
            
            for ii = 1: length(ivec)
                
                if showtext >= 3
                    if round(ii/ceil(length(ivec)/10))*ceil(length(ivec)/10) == ii
                        disp(['      Combination no. ',int2str(ii),' of ',int2str(length(ivec))]) 
                    end
                end
                
                iv = find(jvec==ii);
                ncombs = length(iv);
%                 onesvec3 = ones(ncombs,1);
	
                if showtext >= 4
                    numvec = int2str(edgeshortlist(ii,1));
                    for ll = 2:Ndifforder
                        numvec = [numvec,' ',  int2str(edgeshortlist(ii,ll))];
                    end
                    disp(['   Edges ',numvec])
                        
                end               

                pathalongplane = (edgeseespartialedge(edgeshortlist(ii,2),edgeshortlist(ii,1))>0);
                for jj = 3:Ndifforder
                    pathalongplane = [pathalongplane,(edgeseespartialedge(edgeshortlist(ii,jj),edgeshortlist(ii,jj-1))>0)];                    
                end
                
                if any(multfactors(ivlonglist(iv)))
	
                    IS = full(specextradata(ivlonglist(iv),1:3));
                    IR = full(specextradata(ivlonglist(iv),4:6));
                    
                    firstedgestart = full(edgeextradata(ivlonglist(iv),1));
                    firstedgeend   = full(edgeextradata(ivlonglist(iv),2));
		
                    lastedgestart = full(edgeextradata(ivlonglist(iv),3));
                    lastedgeend   = full(edgeextradata(ivlonglist(iv),4));
		
                    % Calculate the cyl coordinates of all IS/S and IR/R
		
                    cylS = zeros(ncombs,3);
                    edgecoords = [edgestartcoords(edgeshortlist(ii,1),:);edgeendcoords(edgeshortlist(ii,1),:)];
                	[cylS(:,1),cylS(:,2),cylS(:,3)] = EDcoordtrans1(IS,edgecoords,edgenvecs(edgeshortlist(ii,1),:));               
% % % % %                 	[cylS(:,1),cylS(:,2),cylS(:,3)] = EDBcoordtrans1(full(specextradata(ivlonglist(iv),1:3)),edgecoords,edgenvecs(edgeshortlist(ii,1),:));               
		
                    cylR = zeros(ncombs,3);
                    edgecoords = [edgestartcoords(edgeshortlist(ii,Ndifforder),:);edgeendcoords(edgeshortlist(ii,Ndifforder),:)];
                	[cylR(:,1),cylR(:,2),cylR(:,3)] = EDcoordtrans1(IR,edgecoords,edgenvecs(edgeshortlist(ii,Ndifforder),:));     
% % % % %                 	[cylR(:,1),cylR(:,2),cylR(:,3)] = EDBcoordtrans1(full(specextradata(ivlonglist(iv),4:6)),edgecoords,edgenvecs(edgeshortlist(ii,Ndifforder),:));     
                        
                    bc = ones(Ndifforder,2);    % Check real reflfactors!!!    
                    bc = reshape(bc.',2*Ndifforder,1);

                                                
                      if ~isempty(reftoshortlistE)               
                            if Ndifforder >= 2
                                index1 = reftoshortlistE(edgeshortlist(ii,2),edgeshortlist(ii,1));
                                cylE2_r1 = [re1sho(index1,:) thetae1sho(index1,:) ze1sho(index1,:);re2sho(index1,:) thetae2sho(index1,:) ze2sho(index1,:)];
                                index2 = reftoshortlistE(edgeshortlist(ii,1),edgeshortlist(ii,2));
                                cylE1_r2 = [re1sho(index2,:) thetae1sho(index2,:) ze1sho(index2,:);re2sho(index2,:) thetae2sho(index2,:) ze2sho(index2,:)];
                            end
                            if Ndifforder >= 3
                                index1 = reftoshortlistE(edgeshortlist(ii,3),edgeshortlist(ii,2));
                                cylE3_r2 = [re1sho(index1,:) thetae1sho(index1,:) ze1sho(index1,:);re2sho(index1,:) thetae2sho(index1,:) ze2sho(index1,:)];
                                index2 = reftoshortlistE(edgeshortlist(ii,2),edgeshortlist(ii,3));
                                cylE2_r3 = [re1sho(index2,:) thetae1sho(index2,:) ze1sho(index2,:);re2sho(index2,:) thetae2sho(index2,:) ze2sho(index2,:)];
                            end
                            if Ndifforder >= 4
                                index1 = reftoshortlistE(edgeshortlist(ii,4),edgeshortlist(ii,3));
                                cylE4_r3 = [re1sho(index1,:) thetae1sho(index1,:) ze1sho(index1,:);re2sho(index1,:) thetae2sho(index1,:) ze2sho(index1,:)];
                                index2 = reftoshortlistE(edgeshortlist(ii,3),edgeshortlist(ii,4));
                                cylE3_r4 = [re1sho(index2,:) thetae1sho(index2,:) ze1sho(index2,:);re2sho(index2,:) thetae2sho(index2,:) ze2sho(index2,:)];
                            end
                            if Ndifforder >= 5
                                index1 = reftoshortlistE(edgeshortlist(ii,5),edgeshortlist(ii,4));
                                cylE5_r4 = [re1sho(index1,:) thetae1sho(index1,:) ze1sho(index1,:);re2sho(index1,:) thetae2sho(index1,:) ze2sho(index1,:)];
                                index2 = reftoshortlistE(edgeshortlist(ii,4),edgeshortlist(ii,5));
                                cylE4_r5 = [re1sho(index2,:) thetae1sho(index2,:) ze1sho(index2,:);re2sho(index2,:) thetae2sho(index2,:) ze2sho(index2,:)];
                            end
                            if Ndifforder >= 6
                                index1 = reftoshortlistE(edgeshortlist(ii,6),edgeshortlist(ii,5));
                                cylE6_r5 = [re1sho(index1,:) thetae1sho(index1,:) ze1sho(index1,:);re2sho(index1,:) thetae2sho(index1,:) ze2sho(index1,:)];
                                index2 = reftoshortlistE(edgeshortlist(ii,5),edgeshortlist(ii,6));
                                cylE5_r6 = [re1sho(index2,:) thetae1sho(index2,:) ze1sho(index2,:);re2sho(index2,:) thetae2sho(index2,:) ze2sho(index2,:)];
                            end
                            if Ndifforder >= 7
                                index1 = reftoshortlistE(edgeshortlist(ii,7),edgeshortlist(ii,6));
                                cylE7_r6 = [re1sho(index1,:) thetae1sho(index1,:) ze1sho(index1,:);re2sho(index1,:) thetae2sho(index1,:) ze2sho(index1,:)];
                                index2 = reftoshortlistE(edgeshortlist(ii,6),edgeshortlist(ii,7));
                                cylE6_r7 = [re1sho(index2,:) thetae1sho(index2,:) ze1sho(index2,:);re2sho(index2,:) thetae2sho(index2,:) ze2sho(index2,:)];
                            end
                            if Ndifforder >= 8
                                index1 = reftoshortlistE(edgeshortlist(ii,8),edgeshortlist(ii,7));
                                cylE8_r7 = [re1sho(index1,:) thetae1sho(index1,:) ze1sho(index1,:);re2sho(index1,:) thetae2sho(index1,:) ze2sho(index1,:)];
                                index2 = reftoshortlistE(edgeshortlist(ii,7),edgeshortlist(ii,8));
                                cylE7_r8 = [re1sho(index2,:) thetae1sho(index2,:) ze1sho(index2,:);re2sho(index2,:) thetae2sho(index2,:) ze2sho(index2,:)];
                            end
                            
                       else
                            error('STRANGE TO END UP HERE??')
                            if Ndifforder == 2
                                cylE1_r2 = [0 0 edgelengthvec(edgeshortlist(ii,jj));0 0 edgelengthvec(edgeshortlist(ii,jj))]; 
                                cylE2_r1 = [0 0 edgelengthvec(edgeshortlist(ii,jj));0 0 edgelengthvec(edgeshortlist(ii,jj))];
                            else
                                error(['ERROR: Geometries with a single edge can not handle difforder >= 3'])
                            end
                       end

                    for jj = 1:ncombs

                        if multfactors(ivlonglist(iv(jj)))
                            
                            if Ndifforder == 2
                                
                                cylE1_r2frac = cylE1_r2;
                                e1length = edgelengthvec(edgeshortlist(ii,1));
                                cylE2_r1frac = cylE2_r1;
                                e2length = edgelengthvec(edgeshortlist(ii,2));
                                
                                if firstedgestart(jj) ~= 0
                                    cylE1_r2frac(1,3) = cylE1_r2frac(1,3) + e1length*firstedgestart(jj);
                                end
                                if firstedgeend(jj) ~= 1
%                                     cylE1_r2frac(2,3) = cylE1_r2frac(1,3) + e1length*(1-firstedgeend(jj));
                                end
                                if lastedgestart(jj) ~= 0
                                    cylE2_r1frac(1,3) = cylE2_r1frac(1,3) + e2length*lastedgestart(jj);
                                end
                                if lastedgeend(jj) ~= 1
%                                     cylE2_r1frac(2,3) = cylE2_r1frac(1,3)+ e2length*(1-lastedgeend(jj));
                                end
             
                                if midspecpattern(iv(jj))~=0                                    
                                    
                                    % For the cases with specular reflections
                                    % in-between a double diffraction we must
                                    % mirror the two edges to calculate the
                                    % relative-to-edge cylindrical coordinates.
                                    
                                    nspec = sum(midspecpattern(iv(jj),:)>0);
			
                                    % Mirror edge 2 in all the specular reflection
                                    % planes, in reversed order
                                    edgestartmirr = edgestartcoords(edgeshortlist(ii,2),:);
                                    edgeendmirr =   edgeendcoords(edgeshortlist(ii,2),:);
                                    edgevector = edgeendmirr - edgestartmirr;
                                    if lastedgestart(jj) ~= 0
                                        edgestartmirr = edgestartmirr + edgevector*lastedgestart(jj);
                                    else
                                        edgestartmirr = edgestartmirr + edgevector*1e-6;
                                    end
                                    if lastedgeend(jj) ~= 1
                                        edgeendmirr =   edgestartmirr + edgevector*(1-lastedgeend(jj));
                                    else
                                        edgeendmirr =   edgestartmirr + edgevector*(1-1e-6);
                                    end
                                    edgerefcoords = [edgestartcoords(edgeshortlist(ii,1),:);edgeendcoords(edgeshortlist(ii,1),:)];
                                    edgerefnvec =   edgenvecs(edgeshortlist(ii,1),:);
			
                                    % If we have a specular reflection in a plane
                                    % which is perpendicular to the edge plane, we
                                    % should nudge the mirrored edge out a bit so
                                    % that there is no 0/(2*pi) mistake,
                                    if nspec == 1 && ( edgeplaneperptoplane1(midspecpattern(iv(jj),1),edgeshortlist(ii,1)) || edgeplaneperptoplane1(midspecpattern(iv(jj),1),edgeshortlist(ii,2)) )
			%%%%                                reflplanemidpoint = mean(corners(planecorners(midspecpattern(iv(jj),1),1:ncornersperplanevec(midspecpattern(iv(jj),1))),:));
                                        vectowardsmidpoint = approxplanemidpoints(midspecpattern(iv(jj),1),:) - edgestartmirr;
                                        edgestartmirr = edgestartmirr + vectowardsmidpoint*1e-10;
                                        vectowardsmidpoint = approxplanemidpoints(midspecpattern(iv(jj),1),:) - edgeendmirr;
                                        edgeendmirr = edgeendmirr + vectowardsmidpoint*1e-10;
                                    end
                                    xis = [edgestartmirr;edgeendmirr];
                                    for kk = nspec:-1:1
                                        xis = EDBfindis(xis,[midspecpattern(iv(jj),kk);midspecpattern(iv(jj),kk)],planeeqs,2,onesvec1);
                                    end
                                    [rstart,thetastart,zstart,rend,thetaend,zend] = EDBcoordtrans2(xis(1,:),xis(2,:),edgerefcoords,edgerefnvec);
                                    cylE2mirr_r1 = [rstart thetastart zstart;rend thetaend zend];
			
                                    % Mirror edge 1 in all the specular reflection
                                    % planes, in forward order
                                    edgestartmirr = edgestartcoords(edgeshortlist(ii,1),:);
                                    edgeendmirr =   edgeendcoords(edgeshortlist(ii,1),:);
                                    edgevector = edgeendmirr - edgestartmirr;
                                    if firstedgestart(jj) ~= 0
                                        edgestartmirr = edgestartmirr + edgevector*firstedgestart(jj);
                                    else
                                        edgestartmirr = edgestartmirr + edgevector*1e-6;
                                    end
                                    if firstedgeend(jj) ~= 1
                                        edgeendmirr =   edgestartmirr + edgevector*(1-firstedgeend(jj));
                                    else
                                        edgeendmirr =   edgestartmirr + edgevector*(1-1e-6);                               
                                    end
                                    edgerefcoords = [edgestartcoords(edgeshortlist(ii,2),:);edgeendcoords(edgeshortlist(ii,2),:)];
                                    edgerefnvec =   edgenvecs(edgeshortlist(ii,2),:);
			
                                    % Normally, when there is a specular reflection
                                    % in-between, the diffraction path will not be
                                    % along a plane, unless: 
                                    % 1. The same edge is involved twice, with a
                                    % reflection in a perpendicular plane
                                    % in-between.
                                    % 2. See below
                                    pathalongplanewmidspec = 0;
                                    if edgeshortlist(ii,1) == edgeshortlist(ii,2)
                                        if thetastart == 0 || thetastart == (2*pi-closwedangvec(edgeshortlist(ii,2)))
                                            pathalongplanewmidspec = 1;
                                        end
                                    end
                                    
                                    % If we have a specular reflection in a plane
                                    % which is perpendicular to the edge plane, we
                                    % should nudge the mirrored edge out a bit so
                                    % that there is no 0/(2*pi) mistake,
                                    %
                                    % We also have case 2 here (see above) for when
                                    % we could have a pathalongplanewmidspec
                                    % 2. When two different edges have a perpendicular reflection
                                    % in-between 
			
                                    if nspec == 1 && ( edgeplaneperptoplane1(midspecpattern(iv(jj),1),edgeshortlist(ii,1)) || edgeplaneperptoplane1(midspecpattern(iv(jj),1),edgeshortlist(ii,2)) )
                                        % %                                 reflplanemidpoint = mean(corners(planecorners(midspecpattern(iv(jj),1),1:ncornersperplanevec(midspecpattern(iv(jj),1))),:));
                                        vectowardsmidpoint = approxplanemidpoints(midspecpattern(iv(jj),1),:) - edgestartmirr;
                                        edgestartmirr = edgestartmirr + vectowardsmidpoint*1e-10;
                                        vectowardsmidpoint = approxplanemidpoints(midspecpattern(iv(jj),1),:) - edgeendmirr;
                                        edgeendmirr = edgeendmirr + vectowardsmidpoint*1e-10;
                                        pathalongplanewmidspec = 1;
                                    end
                                    xis = [edgestartmirr;edgeendmirr];
                                    for kk = 1:nspec
                                        xis = EDBfindis(xis,[midspecpattern(iv(jj),kk);midspecpattern(iv(jj),kk)],planeeqs,2,onesvec1);
                                    end
                                    [rstart,thetastart,zstart,rend,thetaend,zend] = EDBcoordtrans2(xis(1,:),xis(2,:),edgerefcoords,edgerefnvec);
                                    cylE1mirr_r2 = [rstart thetastart zstart;rend thetaend zend];
                                    
                                    [irnew,~] = EDBwedge2nd(cylS(jj,:),cylR(jj,:),cylE2mirr_r1,cylE1mirr_r2,...
                                        nyvec(edgeshortlist(ii,:)),[edgelengthvec(edgeshortlist(ii,1))*[firstedgestart(jj) firstedgeend(jj)];edgelengthvec(edgeshortlist(ii,2))*[lastedgestart(jj) lastedgeend(jj)]],dzvec(edgeshortlist(ii,:)),...
                                        EDcalcmethod,pathalongplanewmidspec,Rstart,bc,envdata.cair,FSAMP);
                                    
                                else   %  ....   if midspecpattern(iv(jj))~=0,
			
                                    nudge = 0.0000001;

                                    zwedge1 = [edgestartcoords(edgeshortlist(ii,1),:) ; edgeendcoords(edgeshortlist(ii,1),:)];
                                    zwedge2 = [edgestartcoords(edgeshortlist(ii,2),:) ; edgeendcoords(edgeshortlist(ii,2),:)];
                                    z1_end = edgelengthvec(edgeshortlist(ii,1))-nudge;
                                    z2_end = edgelengthvec(edgeshortlist(ii,2))-nudge;

                                    disp('Hello')
                                    savetemp
                                    pause

                                    tfnew = zeros(nfrequencies,1);
                                    for kkk = 1:nfrequencies
%                                         k = 2*pi*frequencies(kkk)/CAIR;
                                        k = 2*pi*frequencies(kkk)/envdata.cair;

                                        tfnew(kkk) =  dblquad('ESIE2betaoverml2_fd',0+nudge,z1_end,0+nudge,z2_end,1e-6,'',k,cylS(jj,:),cylR(jj,:),zwedge1,zwedge2,edgenvecs(edgeshortlist(ii,1),:),edgenvecs(edgeshortlist(ii,2),:),...
                                            nyvec(edgeshortlist(ii,:)),pathalongplane,Rstart,[1 1],envdata.cair);
                                    end
                                        
                                end   %  ....   if midspecpattern(iv(jj))~=0,
			
                                			                
                            elseif Ndifforder >= 3    %   if Ndifforder == 2,
			
                                accuratetripint = 1;
			
                                if accuratetripint == 1   % Sofar this only gives the integral of the irdiff!
                                    
                                    nudge = 1e-7;
			
                                    edge1 = edgeshortlist(ii,1);
                                    edge2 = edgeshortlist(ii,2);
                                    edge3 = edgeshortlist(ii,3);
			
                                    wedgeparams = [cylE2_r1;cylE1_r2;cylE3_r2;cylE2_r3];
                                    edgelengths = [edgelengthvec(edge1) edgelengthvec(edge2) edgelengthvec(edge3)];		
			
                                    tfnew = zeros(nfrequencies,1);
                                    for kkk = 1:nfrequencies
%                                         k = 2*pi*frequencies(kkk)/CAIR;
                                        k = 2*pi*frequencies(kkk)/envdata.cair;
                                        tfnew(kkk) =  triplequad('EDBwedge3rd_intcorefd',0+nudge,edgelengthvec(edge1)-nudge,0+nudge,edgelengthvec(edge2)-nudge,0+nudge,edgelengthvec(edge3)-nudge,1e-9,'',k,cylS(jj,:),cylR(jj,:),wedgeparams,edgelengths,...
                                            edgenvecs(edge1,:),edgenvecs(edge2,:),edgenvecs(edge3,:),...
                                            nyvec(edgeshortlist(ii,:)),pathalongplane,Rstart,[1 1],envdata.cair);
                                    end
                                
                                else  % ...   if accuratetripint == 1,  
			
% % % % %                                     for kk = 1:newndivvec(1)
% % % % %                                         if showtext >= 5
% % % % %                                             disp(['   ',int2str(kk),' of ',int2str(newndivvec(1))]) 
% % % % %                                         end
% % % % %                                         BIGEDGE1stvalue = (kk-0.5)./newndivvec(1);
% % % % %                                         wedgeparams = [cylE2_r1;cylE1_r2;cylE3_r2;cylE2_r3];
% % % % %                                         [irnewpartition,~] = EDBwedgeN(cylS(jj,:),cylR(jj,:),wedgeparams,ncylrows,...
% % % % %                                             nyvec(edgeshortlist(ii,:)),edgelengthvec(edgeshortlist(ii,:)).',...
% % % % %                                             dzvec(edgeshortlist(ii,:)),EDcalcmethod,pathalongplane,nedgeelcombs,Rstart,bc,envdata.cair,FSAMP,BIGEDGE1stvalue);
% % % % %                                         irnewpartition = real(irnewpartition);
% % % % %                                         
% % % % %                                         if kk == 1
% % % % %                                             irnew = irnewpartition;    
% % % % %                                         else
% % % % %                                             lengthaddition = length(irnewpartition);
% % % % %                                             lengthaccum = length(irnew);
% % % % %                                             if lengthaddition > lengthaccum
% % % % %                                                 irnew = [irnew;zeros(lengthaddition-lengthaccum,1)];    
% % % % %                                             end
% % % % %                                             irnew(1:lengthaddition) = irnew(1:lengthaddition) + irnewpartition; 
% % % % %                                         end
% % % % %                                     end
                                    
                                end  % ...   if accuratetripint == 1,  
			
                                if showtext >= 6
                                    disp(['   Calc time = ',num2str(tcalc)])
                                    sum(irnew)
                                    plot(irnew)
                                    pause
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
                            
                            if Ndifforder >= 44   
                                
                                for kk = 1:newndivvec(1)
                                    if showtext >= 5
                                        disp(['   ',int2str(kk),' of ',int2str(newndivvec(1))]) 
                                    end
                                    BIGEDGE1stvalue = (kk-0.5)./newndivvec(1);

                                    [irnewpartition,~] = EDBwedgeN(cylS(jj,:),cylR(jj,:),wedgeparams,ncylrows,...
                                        nyvec(edgeshortlist(ii,:)),edgelengthvec(edgeshortlist(ii,:)).',...
                                        dzvec(edgeshortlist(ii,:)),EDcalcmethod,pathalongplane,nedgeelcombs,Rstart,bc,envdata.cair,FSAMP,BIGEDGE1stvalue);
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
                                                            
                            if showtext >= 5
                                varname = ['allirs',int2str(Ndifforder)];
                                if exist(varname) == 0
                                        eval([varname,' = irnew;'])
                                else
                                        lengthaddition = length(irnew);
                                        eval(['lengthaccum = size(',varname,',1);'])
                                        if lengthaddition > lengthaccum
		                                    eval([varname,' = [',varname,';zeros(lengthaddition-lengthaccum,size(',varname,',2))];'])
                                        end
%                                         allirs2 = [allirs2 zeros(size(allirs2,1),1)];
                                        eval([varname,' = [',varname,' zeros(size(',varname,',1),1)];'])
%                                         allirs2(1:length(irnew),end) = irnew;
                                        eval([varname,'(1:length(irnew),end) = irnew;'])
                                end
                             end
                            
                            % Decide whether the IR should be boosted or not
			
                            boostfactor = 1;
                            if souspecboost ~= 1 || recspecboost ~= 1
                            
                                if prespecpresent(iv(jj)) == 1
                                    boostfactor = souspecboost;    
                                else
                                    if ~isempty(Sinsideplanenumber)
                                        if Sinsideplanenumber(1)~=planesatedge(edgeshortlist(ii,1),1) && Sinsideplanenumber(1)~=planesatedge(edgeshortlist(ii,1),2) 
                                            boostfactor = souspecboost;
                                        end
                                    end    
			
                                end
                                if postspecpresent(iv(jj)) == 1
                                    boostfactor = boostfactor*recspecboost;    
                                else
                                    if ~isempty(Rinsideplanenumber)
                                        if Rinsideplanenumber(1)~=planesatedge(edgeshortlist(ii,Ndifforder),1) && Rinsideplanenumber(1)~=planesatedge(edgeshortlist(ii,Ndifforder),2) 
                                            boostfactor = boostfactor*recspecboost;
                                        end
                                    end
                                end
                            
                            end   % ...   if souspecboost ~= 1 | recspecboost ~= 1,
                            
                           % 21 Oct. 2017
                            % The claim on the three lines below is not
                            % always correct!
                                % For thin plates, we must have a boost factor!
                                % This is because there will be multiple equivalent
                                % combinations passing on the rear side of the thin plate
                            
                                % Old code
%                             if all( nyvec(edgeshortlist(ii,:)) == 0.5 )
%                                 boostfactor = boostfactor*2^(Ndifforder-1);                        
%                             end

                            % New code 21 Oct. 2017
                            
                            if length(pathalongplane) == 1
                                if pathalongplane == 1
                                    boostfactor = boostfactor*2;
                                end                                
                            else
                               error('ERROR: Correct this: inrease boost factor for all paths that travel along a thin plate'); 
                            end
                            
                            if boostfactor ~= 1
                                tfnew = tfnew*boostfactor;    
                            end
                            
                            tfnew = tfnew*double(multfactors(ivlonglist(iv(jj))));
			
%                             if saveindividualdifftfs(1) == 0,
                                tfdiff = tfdiff + tfnew;
%                             else
%                                 ncolsdiff = size(tfdiff,2);
%                                 if Ndifforder > ncolsdiff,
%                                     tfdiff = [tfdiff zeros(nfrequencies,Ndifforder-ncolsdiff)];    
%                                     ncolsdiff = ncolsdiff+1;
%                                 end
%                                 tfdiff(:,Ndifforder) = tfdiff(:,Ndifforder) + tfnew;
%                             end
                            
                        end    % ...  if multfactors(ivlonglist(iv(jj))),
                        
                    end    % ....for jj = 1:ncombs,

                else   %  ...    if any(multfactors(ivlonglist(iv))),
                    if showtext >= 4
                        disp('      Combination not computed because of repetition')                         
                    end
                    
                end    % ...    if any(multfactors(ivlonglist(iv))),
                
            end   % .... for ii = 1: length(ivec),
        
        end   % ...     if any(difforderinlist==Ndifforder),
    
%         if size(tfdiff,2) < Ndifforder & saveindividualdifftfs(1) == 1,
%             tfdiff = [tfdiff zeros(nfrequencies,Ndifforder-size(tfdiff,2))];    
%         end
	end    % ... for Ndifforder = 2:specorder,
end     % .... if userwantsdiff2 == 1,
