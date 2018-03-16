function [irhod,EDinputdatahash] = EDmakeHODirs(hodpaths,difforder,elemsize,edgedata,...
    edgetoedgedata,Sdata,doaddsources,sourceamplitudes,Rdata,cair,fs,Rstart,showtext,EDversionnumber)
% EDmakeHODirs - Constructs higher-order (two and higher) diffraction impulse
% responses from a list of paths in the input struct hodpaths.
%
% Input parameters:
%   hodpaths        Cell variable with a list of all the valid edge-combos
%                   for each diffraction order
%   difforder       Value from controlparameters
%   elemsize        Vector, [1,difforder], with values that specify how
%                   fine edge discretization will be used for each order.
%                   Recommended starting values are [1 0.5 0.25 0.125 ...].
%                   Higher values give more accurate results.
%   edgedata,edgetoedgedata,Sdata,Rdata     Structs
%   doaddsources    0 or 1; determines if individual irs for each source
%                   should be stored, or if they should be summed
%   sourceamplitudes   the amplitude that will be used when all sources' 
%                   contributions are added
%   cair,fs,Rstart  Values from envdata and controlparameters 
%   showtext        Value from filehandlingparameters
%   EDversionnumber
%
% Output parameters:
%   irhod           The ir with all the higher-order diffraction orders summed up
%                   Size:
%                       [nsampels,nreceivers,nsources] (if doaddsources = 0)
%                       [nsampels,nreceivers,1] (if doaddsources = 1)
%   EDinputdatahash
%
% Uses functions EDcreindexmatrix,EDwedge2nd,EDwedgeN from the EDtoolbox
% Uses function DataHash from Matlab Central
%
% Peter Svensson (peter.svensson@ntnu.no) 16 Mar 2018
%
% [irhod,EDinputdatahash] = EDmakeHODirs(hodpaths,difforder,elemsize,edgedata,...
%     edgetoedgedata,Sdata,doaddsources,sourceamplitudes,Rdata,cair,fs,Rstart,showtext);
    
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

global BIGEDGESTEPMATRIX 

EDinputdatastruct = struct('difforder',difforder,'hodpaths',hodpaths,'elemsize',elemsize,...
'edgedata',edgedata,'Sdata',Sdata,'doaddsources',doaddsources,'sourceamplitudes',sourceamplitudes,'Rdata',Rdata,'cair',cair,...
'fs',fs,'Rstart',Rstart,'EDversionnumber',EDversionnumber);
EDinputdatahash = DataHash(EDinputdatastruct);

%--------------------------------------------------------------------------

nyvec = pi./(2*pi - edgedata.closwedangvec);    

reftoshortlistE = edgetoedgedata.reftoshortlistE;

nsources = size(Sdata.sources,1);
nreceivers = size(Rdata.receivers,1);

firstcomponentdone = 0;

%--------------------------------------------------------------------------

for isou = 1:nsources
    for irec = 1:nreceivers

        S = Sdata.sources(isou,:);
        R = Rdata.receivers(irec,:);

        for Ndifforder = 2:difforder
            if showtext >= 2
                disp(['   Diffraction order ',int2str(Ndifforder)])    
            end

            if ~isempty(hodpaths{Ndifforder})

                pathalongplane = ones(1,Ndifforder-1);
                bc = ones(1,Ndifforder);

                % Calculate some general parameters that are shared for all
                % N-diffraction calculations

                divmin = cair/(fs*elemsize(Ndifforder));
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

                    if showtext >= 3
                        if round(ii/ceil(ncomponents/10))*ceil(ncomponents/10) == ii
                            disp(['      Combination no. ',int2str(ii),' of ',int2str(ncomponents)]) 
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
                            'n',pathalongplane,Rstart,bc,cair,fs);

        %                 IRDIFFVEC = [IRDIFFVEC;sum(irnew)];

                    elseif Ndifforder == 3    %   if Ndifforder == 2			

                        for kk = 1:newndivvec(1)                    
                            BIGEDGE1stvalue = (kk-0.5)./newndivvec(1);
                            wedgeparams = [cylE2_r1;cylE1_r2;cylE3_r2;cylE2_r3];

                            [irnewpartition,~] = EDwedgeN(cylS,cylR,wedgeparams,ncylrows,...
                                nyvec(edgepattern(:)),edgedata.edgelengthvec(edgepattern(:)).',...
                                dzvec(edgepattern(:)),'n',pathalongplane,nedgeelcombs,Rstart,bc,cair,fs,BIGEDGE1stvalue);
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
                                dzvec(edgepattern(:)),'n',pathalongplane,nedgeelcombs,Rstart,bc,cair,fs,BIGEDGE1stvalue);
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

            %                             if all( nyvec(edgepattern(:)) == 0.5 )
                   boostfactor = 1;
                   if all( nyvec(edgepattern(:)) == 0.5 ) && pathalongplane                                
                        boostfactor = boostfactor*2^(Ndifforder-1);                        
                    end

                    if boostfactor ~= 1
                        irnew = irnew*boostfactor;    
                    end

                    nnew = length(irnew);

                    if firstcomponentdone == 0
                       if doaddsources == 0
                           irhod = zeros(nnew,nreceivers,nsources);
                           irhod(:,irec,isou) = irnew;
                       else
                           irhod = zeros(nnew,nreceivers,1);
                           irhod(:,irec,1) = irnew*sourceamplitudes(isou);                   
                       end
                       firstcomponentdone = 1;
                    else
                        nold = size(irhod,1);
                        if nnew > nold
                            if doaddsources == 0
                               irhod = [irhod;zeros(nnew-nold,nreceivers,nsources)]; 
                            else
                               irhod = [irhod;zeros(nnew-nold,nreceivers,1)];                         
                            end
                        end
                        if doaddsources == 0
                            irhod(1:nnew,irec,isou) = irhod(1:nnew,irec,isou) + irnew;                
                        else
                            irhod(1:nnew,irec,1) = irhod(1:nnew,irec,1) + irnew*sourceamplitudes(isou);                                    
                        end
                    end

                end 

            end  

        end  

    end
    
end
