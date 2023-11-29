function[tfinteqdiff,timingdata,extraoutputdata,elapsedtimehodtf,existingfilename] = ...
    EDintegralequation_convex_tf(Hsubmatrixdata,planedata,edgedata,edgetoedgedata,...
	Sdata,Rdata,envdata,controlparameters,EDversionnumber,filehandlingparameters)
% EDintegralequation_convex_tf calculates the sound pressure representing second-
% and higher-order diffraction, for a convex scattering object
%
% From v0.4 of the EDtoolbox, the input parameter list changed.
% 
% Input parameters:
%   Hsubmatrixdata          Structs
%   planedata,edgedata,edgetoedgedata       Structs
%   Sdata,Rdata,envdata,controlparameters   Structs
%   EDversionnumber
%   filehandlingparameters, envdata, 
%
% Output parameters
%   tfinteqdiff      Matrix with the diffracted pressure
%                   If nsources = 1, or all sources should be added, then
%                       the size is [nfrequencies,nreceivers]
%                   If nsources > 1, and all sources should be computed
%                       separately, then the size is
%                       [nfrequencies,nreceivers,nsources]
%   timingdata
%   extraoutputdata
%   elapsedtimehodtf
%                       This tells how long time was used inside this
%                       function. If an existing file was reused, then
%                       elapsedtimehodtf has a second value which tells
%                       how much time was used for the existing file.
%   existingfilename    For v2 of this function, if an existing file was
%                       found that was reused, then the reused file
%                       name is given here. If no existing file could be 
%                       reused then this variable is empty. For v1 of
%                       this function, this variable is also empty.
%
% Uses functions EDdistelements, EDcalcedgeinteqmatrixsub2_mex, EDinteg_souterm, 
% EDcalcpropagatematrix, EDrecycleresultfiles, EDcoordtrans1 from EDtoolbox
% Uses function DataHash from Matlab Central
%           
% Peter Svensson (peter.svensson@ntnu.no)  29 Nov. 2023  
%                       
% [tfinteqdiff,timingdata,extraoutputdata,elapsedtimehodtf,existingfilename] = ...
%     EDintegralequation_convex_tf(Hsubmatrixdata,planedata,edgedata,edgetoedgedata,
%	Snewdata,Rnewdata,envdata,controlparameters,EDversionnumber,filehandlingparameters);

% 31 March 2015 Introduced detailed timing, also as output parameter.
% 8 April 2015 Substantial speeding up by saving Hsubdata instead of Hsub
%              and using the accumarray function.
% 28 June 2015 Added the option of saving the edge source signals.
% 29 Jan. 2016 Introduced the input parameters calcinteq_souterms,calcinteq_edgeterms,calcinteq_propagate
% 30 Jan 2016 Changed the inteqsettings to a struct
% 29 April 2016 Fixed a small bug (niterations did not use the struct
%              correctly in a showtext case). Also corrected a bug with
%              edge-to-edge transfer between non-parallel edges.
% 2 March 2017 Introduced a swtich between mac and non-mac computers, since
%              a mex-file exists only for mac.
%              Also, corrected a use of uint8 such that more than 255 edge
%              points are possible.
% 10 March 2017 Removed the automatic saving of ivproblematic; now it will
%               happen only if ivproblematicmatrix contains non-zero values.
% 10 March 2017 Introduced doaddsources = 2 for regular source/receiver
%               matrices which exploits symmetries.
% 25 April 2017 A small change in the for loop over the Hsub matrices. A
%               list is constructed so that the iterative loop doesnt have
%               to check if the Hsub exists or not.
% 30 May 2017   Changed back what was changed on 25 April, because it was
%               wrong.
% 26 Oct. 2017 Fixed a bug that occurred for skewed edges: the r-values
%              were distributed linearly before.
% 4 Nov. 2017  Implemented the partly visible edge in the source term.
% 27 Nov. 2017 Copied from ESIE2toolbox and trimmed down + changed to the
%               use of cells for the Hsubs.
% 28 Nov. 2017 Cleaned up code a bit. Added timingdata.
% 29 Nov. 2017 Cleaned up: removed the old timingdata code. Reduced the
%              no. of displayed frequencies.
% 13 Dec. 2017 Added the sourceamplitudes to the Sinputdata struct
% 19 Jan 2018 Renamed the ESIE2...mex to ED...mex

% 10 Dec. 2017 TODO: finish the timingdata for the CASE2 : addsources
% 19 Jan 2018 The TODO above is probably done
% 23 Jan 2018 Changed to "~strcmp(computer, 'MACI64')" after hint from Jan
%             Slechta.
% 26 Jan 2018 Implemented the handling of doallSRcombinations = 0.
% 31 Jan 2018 Had not implemented sourceamplitudes for case 2.
% 9 Feb 2018 Introduced EDinputdatahash
% 28 Feb 2018 The saving of sousig files had an extra "results"
% subdirectory.
% 13 Apr 2018 Fixed a small bug on l143: thinplaneboostvec got the format
% uint8 but needs to be double. Found by Antoine.
% 3 June 2020 Fixed a bug: folder names with spaces can be handled now
% 28 Sep. 2023 Implemented version 2 of this function while maintaining
% compatibility with the old "version 1". v2 moves the check if an existing
% file can be reused inside this function. Also updated load and save to
% the function call form, which avoids problems with spaces in file names.
% 3 Oct. 2023 Introduced the EDsettingshash as input parameter.
% 6 Oct. 2023 Fixed a small error: the timing data wasn't returned when an
% existing file was reused.
% 27 Oct. 2023 Changed the input parameter list substantially.
% 30 Oct. 2023 Fine-tuned the EDinputdatahash
% 1 Nov. 2023 Fixed a small error: doaddsources and sourceamplitudes and 
% doallSRcombinations were missing the Sdata. in one location.
% 20 Nov. 2023 Implemented the piston source
% 24 Nov. 2023 Stores many gaussvalues, up to the max. needed
% 29 Nov. 2023 Adapted to the name change of the field pistongausspoints to
% pistongaussorder

t00 = clock;
showtext = filehandlingparameters.showtext;

EDinputdatastruct = struct('corners',planedata.corners,'planecorners',...
    planedata.planecorners,'offedges',edgedata.offedges,...
    'sources',Sdata.coordinates,'Snedgesubs',Sdata.nedgesubs,...
    'receivers',Rdata.coordinates,'Rnedgesubs',Rdata.nedgesubs,...
    'difforder',controlparameters.difforder,'frequencies',...
    controlparameters.frequencies,'ngauss',controlparameters.ngauss, ...
    'discretizationtype',controlparameters.discretizationtype,...
    'Rstart',controlparameters.Rstart,'EDversionnumber',EDversionnumber);
EDinputdatahash = DataHash(EDinputdatastruct);

%---------------------------------------------------------------
% Sort out the file business: can an existing file be used?
% Then copy the existing file to a new copy. 

if filehandlingparameters.suppressresultrecycling == 1
    foundmatch = 0;
    existingfilename = '';
else
    [foundmatch,existingfilename] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_tfinteq',EDinputdatahash);
end

EDinputdatahash_tfinteq = EDinputdatahash;

desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tfinteq.mat'];

if foundmatch == 1
    eval(['load(''',existingfilename,''')'])
    if ~strcmp(existingfilename,desiredname)
        copyfile(existingfilename,desiredname);
    end
    elapsedtimehodtf_new = etime(clock,t00);
    elapsedtimehodtf = [elapsedtimehodtf_new elapsedtimehodtf];
    P_receiver = tfinteqdiff;
    return
end

extraoutputdata = struct('tfinteqdiff_nodiff2',[]);

nfrequencies = length(controlparameters.frequencies);

timingdata = zeros(1,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detect special cases: check if the object is thin
    
nsources = size(Sdata.coordinates,1);
nreceivers = size(Rdata.coordinates,1);

% Find if there are image sources that could see any edges
% (NB! This is not relevant for convex scattering objects)
nplanes = size(planedata.planecorners,1);
nedges = size(edgedata.edgecorners,1);
vispartedgesfromIS = sparse(zeros(nplanes,nedges));

gaussvectors = cell(controlparameters.ngauss,1);
for ii = 2:controlparameters.ngauss
    [gaussvec1,gaussvec2] = EDdistelements(ii,controlparameters.discretizationtype);
    gaussvectors{ii} = [gaussvec1 gaussvec2];
end

ivproblematicmatrix = zeros(nreceivers,100);

listofactiveHsubs = zeros(size(Hsubmatrixdata.listofsubmatrices,1),1);
Hsubcounter = 1;

alltripletsthin = prod(double(Hsubmatrixdata.isthinplanetriplet));
notripletsthin = prod(1-double(Hsubmatrixdata.isthinplanetriplet));
nbig = Hsubmatrixdata.bigmatrixendnums(end);
nsubmatrices = size(Hsubmatrixdata.listofsubmatrices,1);

refto = edgetoedgedata.reftoshortlistE;

if isempty(Hsubmatrixdata.reftoshortlist)
    symmetrycompression = 0;
else    
    symmetrycompression = 1;
end

if symmetrycompression < 1
   reftoshortlist = 1:size(Hsubmatrixdata.listofsubmatrices,1);
end

Hsubconstant = -1/4/pi/2;

if alltripletsthin == 1
    thinplaneboostvec = 2*ones(nsubmatrices,1);
else
    if notripletsthin == 1
        thinplaneboostvec = ones(nsubmatrices,1);
    else
        thinplaneboostvec = double(Hsubmatrixdata.isthinplanetriplet(1:nsubmatrices)+1);
    end
end

Hsubdatalists = cell(nsubmatrices,1);
ivuselists = cell(nsubmatrices,1);

isthinhole = 0;

% ----------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------
% We must have an outer loop over frequency, since all matrices change with
% frequency. The for-loop stretches to the end of this function.

if showtext >= 1
    freqstep = round(nfrequencies/10);
    if freqstep == 0
       freqstep = 1; 
    end
end

for ifreq = 1:nfrequencies     
    frequency = controlparameters.frequencies(ifreq);
    if showtext >= 1
        if round(ifreq/freqstep)*freqstep == ifreq
            disp(['      Frequency no. ',int2str(ifreq),' (of ',int2str(nfrequencies),'): ',num2str(frequency),' Hz'])
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PART 1: Build the H-matrix

    if showtext >= 2 && controlparameters.difforder > 2
        disp(' ')
        disp(['      Building the H matrix of size (',int2str(nbig),',',int2str(nbig),') with ',int2str(Hsubmatrixdata.submatrixcounter),' submatrices'])
    end
    
    if ifreq == 1
        t00 = clock;
    end

    k = 2*pi*frequency/envdata.cair;
    
    % In the for-loop, nn will be the Hsub number
    % The ii value will refer to the row-number in the edgepairlist for the 
    % receiving edge-pair.
    % The jj value will refer to the row-number in the edgepairlist for the 
    % sending edge-pair.

    n1prev = 0;
    n2prev = 0;
    n3prev = 0;
    shortlistprev = 0;

    % ----------------------------------------------------------------------------------------------
    % ----------------------------------------------------------------------------------------------
    % ----------------------------------------------------------------------------------------------
    % If difforder > 2, create the big H-matrix, in the form of separate sub-matrices
    % The H-matrix is independent of sources and receivers.

    if controlparameters.difforder > 2
        for nn = 1:nsubmatrices

            edge1 = Hsubmatrixdata.edgetripletlist(nn,3);
            n1 = Hsubmatrixdata.nedgeelems(edge1);

            edge2 = Hsubmatrixdata.edgetripletlist(nn,2);                   
            n2 = Hsubmatrixdata.nedgeelems(edge2);
               len2 = edgedata.edgelengthvec(edge2);

            edge3 = Hsubmatrixdata.edgetripletlist(nn,1);
            n3 = Hsubmatrixdata.nedgeelems(edge3);
%             disp(['triplet no. ',int2str(nn),': n1,n2,n3 = ',int2str(n1),',',int2str(n2),',',int2str(n3)])
            shortlistnumber = Hsubmatrixdata.reftoshortlist(nn);

%             ind_start = 1 + (Hsubmatrixdata.bigmatrixstartnums(ii)-1) + (Hsubmatrixdata.bigmatrixstartnums(jj)-1)*nbig;

           if showtext >= 3
               disp(' ')
                disp(['         H-submatrix ',int2str(nn),': From edge ',int2str(edge1),' via edge ',int2str(edge2),' to edge ',int2str(edge3)])
                disp(['         of size ',int2str(n3*n2),' by ',int2str(n2*n1),'. Thinplaneboost = ',int2str(thinplaneboostvec(nn))])
                if symmetrycompression ==1
                    disp(['         stored as shortlist number ',int2str(shortlistnumber)])
                end
           end

           % Build vertical [n3,n2] matrices and expand them horizontally
           % NB! n1 refers to edge 1, which is the source edge
           %     n3 refers to edge 3, which is the target edge
           % edgetriplets written in order [edge3 edge2 edge 1]

           if n1 ~= n1prev || n2 ~= n2prev || n3 ~= n3prev

%                if computer ~= 'MACI64'
               if ~strcmp(computer, 'MACI64')
                   if n3 <= 255
                       n3vertmat = uint8(1:n3);
                   else
                       n3vertmat = uint16(1:n3);                   
                   end
                   n3vertmat = reshape(repmat(n3vertmat,n2,1),n2*n3,1);           
                   n3vertmat = n3vertmat(:,ones(1,n1*n2));

                   if n2 <= 255
                       n2vertmat = uint8(1:n2).';
                   else
                      n2vertmat = uint16(1:n2).'; 
                   end
                   n2vertmat = repmat(n2vertmat,n3,1);
                   n2vertmat = n2vertmat(:,ones(1,n1*n2));

               % Build horizontal n2 and n1 matrices and expand them vertically
                    if n1 <= 255                    
                       n1hormat = uint8(1:n1);
                    else
                        n1hormat = uint16(1:n1);
                    end
                   n1hormat = repmat(n1hormat,1,n2);
                       n1hormat = n1hormat(ones(n3*n2,1),:);

                    if n2 <= 255
                       n2hormat = uint8(1:n2);
                    else
                       n2hormat = uint16(1:n2); 
                    end
                       n2hormat = reshape(repmat(n2hormat,n1,1),1,n2*n1);           
                   n2hormat = n2hormat(ones(n3*n2,1),:);

               % Now, we should use only the combinations where n2vertmat ==
               % n2hormat! Otherwise we get too many combinations: there should
               % be n1*n2*n3 unique combinations, but the Hsub-matrix has size
               % [n3*n2,n2*n1], which has n2 too many elements!
               %
               % Once we have the iv, which identifies the wanted matrix
               % elements, we can throw away n2hormat because it is identical
               % to n2vertmat in the wanted matrix elements!

                   ivuse = find(n2hormat==n2vertmat);
                   n3vertmat = n3vertmat(ivuse);
                   n2vertmat = n2vertmat(ivuse);
                   n1hormat  = n1hormat(ivuse);

               else
%                     [ivuse,n1hormat,n2vertmat,n3vertmat] = ESIE2calcedgeinteqmatrixsub2_mex(n1,n2,n3);       
                    [ivuse,n1hormat,n2vertmat,n3vertmat] = EDcalcedgeinteqmatrixsub2_mex(n1,n2,n3);       

               end

               ivuse = uint64(ivuse);

           end

           % Build vectors that will be referred to by the n1,n2,n3
           % matrices

           if n1 ~= n1prev
                n1vec = full(Hsubmatrixdata.quadraturematrix_pos(n1,1:n1));
                n1vec = n1vec(:);
                weightvec = full(Hsubmatrixdata.quadraturematrix_weights(n1,1:n1));
                weightvec = weightvec(:);
                n1prev = n1;
           end       
           dzvec = weightvec*edgedata.edgelengthvec(edge1);

            if n2 ~= n2prev
                n2vec = full(Hsubmatrixdata.quadraturematrix_pos(n2,1:n2));
                n2vec = n2vec(:);
                n2prev = n2;
            end
            if n3 ~= n3prev        
                n3vec = full(Hsubmatrixdata.quadraturematrix_pos(n3,1:n3));
                n3vec = n3vec(:);
                n3prev = n3;
            end

            % Bug found and fixed 22 Sep. 2014: lines below did not always detect if another
            % edge was skewed relative to the first one, because only the
            % thetae1sho was used. Necessary to check thetae1sho and thetae2sho.

           ny = pi/(2*pi-edgedata.closwedangvec(edge2));
           thetain  = edgetoedgedata.thetae1sho(refto(edge1,edge2));
           thetain_other_end_of_edge  = edgetoedgedata.thetae2sho(refto(edge1,edge2));
        %    if abs(thetain) < 1e-10 || abs(thetain-(2*pi-closwedangvec(edge2))) < 1e-10   % Old erroneous line
           if (abs(thetain) < 1e-10 || abs(thetain-(2*pi-edgedata.closwedangvec(edge2))) < 1e-10) && abs(thetain-thetain_other_end_of_edge) < 1e-10
               acrossface_in = 1;       
           else
               acrossface_in = 0;
           end
           thetaout = edgetoedgedata.thetae1sho(refto(edge3,edge2));
           thetaout_other_end_of_edge = edgetoedgedata.thetae2sho(refto(edge3,edge2));
        %   if abs(thetaout) < 1e-10 || abs(thetaout-(2*pi-closwedangvec(edge2))) < 1e-10    % Old erroneous line
           if (abs(thetaout) < 1e-10 || abs(thetaout-(2*pi-edgedata.closwedangvec(edge2))) < 1e-10) && abs(thetaout-thetaout_other_end_of_edge) < 1e-10
               acrossface_out = 1;
           else
               acrossface_out = 0;
           end

           if showtext >= 3
               disp(['         Acrossface-out = ',int2str(acrossface_out),' and Acrossface-in = ',int2str(acrossface_in), ' (not used)'])
           end

           ze2_re2 = n2vec*len2;
           ze3_re2 = edgetoedgedata.ze1sho(refto(edge3,edge2)) + n3vec*( edgetoedgedata.ze2sho(refto(edge3,edge2))-edgetoedgedata.ze1sho(refto(edge3,edge2))   );
    %        re3_re2 = edgetoedgedata.re1sho(refto(edge3,edge2)) + n3vec*( edgetoedgedata.re2sho(refto(edge3,edge2))-edgetoedgedata.re1sho(refto(edge3,edge2))   );

           if shortlistnumber ~= shortlistprev      

               n1vec = n1vec(:).';

               ze1_re2 = edgetoedgedata.ze1sho(refto(edge1,edge2)) + n1vec*( edgetoedgedata.ze2sho(refto(edge1,edge2))-edgetoedgedata.ze1sho(refto(edge1,edge2))   );

               % Bug fixed 26 Oct. 2017: the line below is wrong because the
               % r-values are not distributed linearly.
    %            re1_re2 = edgetoedgedata.re1sho(refto(edge1,edge2)) + n1vec*( edgetoedgedata.re2sho(refto(edge1,edge2))-edgetoedgedata.re1sho(refto(edge1,edge2))   );


               if acrossface_in*acrossface_out == 0 % This case is a preparation for non-convex objects.
                    % Reconvert the cylindrical coordinates of
                    % edge1-re2 into cartesian! Then we can use the EDcoordtrans.
                    % We define our own cartesian coord syst such that the reference
                    % edge (edge2) has its starting point in [0 0 0], the z-axis along its
                    % own edge, and the x-axis is along (onto) the reference plane.

                    xedge1start = edgetoedgedata.re1sho(refto(edge1,edge2))*cos(edgetoedgedata.thetae1sho(refto(edge1,edge2)));
                    xedge1end   = edgetoedgedata.re2sho(refto(edge1,edge2))*cos(edgetoedgedata.thetae2sho(refto(edge1,edge2)));
                    yedge1start = edgetoedgedata.re1sho(refto(edge1,edge2))*sin(edgetoedgedata.thetae1sho(refto(edge1,edge2)));
                    yedge1end   = edgetoedgedata.re2sho(refto(edge1,edge2))*sin(edgetoedgedata.thetae2sho(refto(edge1,edge2)));
                    xe1_re2 = xedge1start + n1vec.*(xedge1end-xedge1start);
                    ye1_re2 = yedge1start + n1vec.*(yedge1end-yedge1start);

                    [re1_re2,thetae1_re2,~] = EDcoordtrans1([xe1_re2.' ye1_re2.' ze1_re2.'],[0 0 0;0 0 edgedata.edgelengthvec(edge2)],[0 1 0]);
                    re1_re2 = re1_re2(:).';
                    thetae1_re2 = thetae1_re2(:).';

                    xedge3start = edgetoedgedata.re1sho(refto(edge3,edge2))*cos(edgetoedgedata.thetae1sho(refto(edge3,edge2)));
                    xedge3end   = edgetoedgedata.re2sho(refto(edge3,edge2))*cos(edgetoedgedata.thetae2sho(refto(edge3,edge2)));
                    yedge3start = edgetoedgedata.re1sho(refto(edge3,edge2))*sin(edgetoedgedata.thetae1sho(refto(edge3,edge2)));
                    yedge3end   = edgetoedgedata.re2sho(refto(edge3,edge2))*sin(edgetoedgedata.thetae2sho(refto(edge3,edge2)));
                    xe3_re2 = xedge3start + n3vec.*(xedge3end-xedge3start);
                    ye3_re2 = yedge3start + n3vec.*(yedge3end-yedge3start);

                    [re3_re2,thetae3_re2,~] = EDcoordtrans1([xe3_re2 ye3_re2 ze3_re2],[0 0 0;0 0 edgedata.edgelengthvec(edge2)],[0 1 0]);
                    re3_re2 = re3_re2(:);
                    thetae3_re2 = thetae3_re2(:);

               else

                    re1_re2 = edgetoedgedata.re1sho(refto(edge1,edge2)) + n1vec*( edgetoedgedata.re2sho(refto(edge1,edge2))-edgetoedgedata.re1sho(refto(edge1,edge2))   );
                    re3_re2 = edgetoedgedata.re1sho(refto(edge3,edge2)) + n3vec*( edgetoedgedata.re2sho(refto(edge3,edge2))-edgetoedgedata.re1sho(refto(edge3,edge2))   );

               end

               ldist = (ze3_re2(n3vertmat) - ze2_re2(n2vertmat)).^2   + re3_re2(n3vertmat).^2;
               ldist = sqrt(ldist);

               mdist = (ze1_re2(n1hormat).' - ze2_re2(n2vertmat)).^2   + (re1_re2(n1hormat).').^2;
               mdist = sqrt(mdist);

               dz_expjkmoverm = exp(-1i*k*mdist)./mdist.*dzvec(n1hormat);           

               ch = ((ze2_re2(n2vertmat) - ze1_re2(n1hormat).').*(ze2_re2(n2vertmat) - ze3_re2(n3vertmat)) + mdist.*ldist)./re1_re2(n1hormat).'./re3_re2(n3vertmat);
               ch = ( real(sqrt( ch.^2-1)) + ch ).^ny;
               ch = ( ch + 1./ch)/2;

               % PS 28 March 2011: The 4 beta terms are identical!       
               if acrossface_in == 1 && acrossface_out == 1
                   beta = 4*( sin(ny*(pi + thetain + thetaout  ))./(ch - cos(ny*(pi + thetain + thetaout   ))));
               else
                   beta = ( sin(ny*(pi + thetae1_re2(n1hormat).' + thetae3_re2(n3vertmat)  ))./(ch - cos(ny*(pi + thetae1_re2(n1hormat).' + thetae3_re2(n3vertmat)   ))+eps*10)) + ...
                          ( sin(ny*(pi + thetae1_re2(n1hormat).' - thetae3_re2(n3vertmat)  ))./(ch - cos(ny*(pi + thetae1_re2(n1hormat).' - thetae3_re2(n3vertmat)   ))+eps*10)) + ...
                          ( sin(ny*(pi - thetae1_re2(n1hormat).' + thetae3_re2(n3vertmat)  ))./(ch - cos(ny*(pi - thetae1_re2(n1hormat).' + thetae3_re2(n3vertmat)   ))+eps*10)) + ...
                          ( sin(ny*(pi - thetae1_re2(n1hormat).' - thetae3_re2(n3vertmat)  ))./(ch - cos(ny*(pi - thetae1_re2(n1hormat).' - thetae3_re2(n3vertmat)   ))+eps*10));
               end

               if ismember(edge1,edgedata.offedges)==1 || ismember(edge2,edgedata.offedges)==1 || ismember(edge3,edgedata.offedges)==1
                   beta = beta*0;
               else
        %           disp(['Edge1: ',int2str(edge1),' edge2: ',int2str(edge2),' edge3: ',int2str(edge3)])
               end

                Hsubdatalists{nn} = (2-acrossface_out)*Hsubconstant*ny*beta.*dz_expjkmoverm*thinplaneboostvec(nn);   
                ivuselists{nn} = ivuse;
                if ifreq == 1
                    listofactiveHsubs(Hsubcounter) = nn;
                    Hsubcounter = Hsubcounter + 1;
                end

           end

            shortlistprev = shortlistnumber;

        end
    end
    
    if ifreq == 1
        ivactive = listofactiveHsubs~=0;
        listofactiveHsubs = listofactiveHsubs(ivactive);
    end
    
    if ifreq == 1
        timingdata(1) = etime(clock,t00);
        t00 = clock;
    end

    % End of PART 1: All the H submatrices have been built (for one freq.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CASE 1: A single source, or many sources' contributions should be
    % added to give the the edge source signals
    % (doaddsources = 1, or we have only one source)
    % PART 2: Create Q_firstterm 
    
    if Sdata.doaddsources == 1 || nsources == 1
        if ifreq==1
            P_receiver = zeros(nfrequencies,nreceivers);
        end
        
        if filehandlingparameters.loadinteqsousigs == 0
                        
            Q_firstterm = 0;
            for isou = 1:nsources
                ISOU = int2str(isou);

                if strcmp(Sdata.sourcetype,'monopole') == 1
                     Q_firstterm_addition = EDinteg_souterm(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,...                    
                        controlparameters,Sdata.vispartedgesfroms(:,isou),Sdata.vispartedgesfroms_start(:,isou),...
                        Sdata.vispartedgesfroms_end(:,isou),frequency,gaussvectors,Sdata.rSsho(Sdata.reftoshortlistS(:,isou)),...
                        Sdata.thetaSsho(Sdata.reftoshortlistS(:,isou)),Sdata.zSsho(Sdata.reftoshortlistS(:,isou)),filehandlingparameters.showtext);
                elseif strcmp(Sdata.sourcetype,'polygonpiston') == 1
                    Q_firstterm_sumpistongaussorder = 0;
                    for np = 1:size(Sdata.pistongausscoordinates{1},1)
                        Q_firstterm_addition = EDinteg_souterm(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,...                    
                            controlparameters,Sdata.vispartedgesfroms(:,isou),Sdata.vispartedgesfroms_start(:,isou),...
                            Sdata.vispartedgesfroms_end(:,isou),frequency,gaussvectors,Sdata.rSsho(Sdata.reftoshortlistS(:,isou,np)),...
                            Sdata.thetaSsho(Sdata.reftoshortlistS(:,isou,np)),Sdata.zSsho(Sdata.reftoshortlistS(:,isou,np)),filehandlingparameters.showtext);
                        Q_firstterm_sumpistongaussorder = Q_firstterm_sumpistongaussorder + Q_firstterm_addition*squeeze(Sdata.pistongaussweights{isou}(np,1));
                    end
                    Q_firstterm_addition = Q_firstterm_sumpistongaussorder;
                end
%                 Q_firstterm = Q_firstterm + Q_firstterm_addition;
                Q_firstterm = Q_firstterm + Q_firstterm_addition*Sdata.sourceamplitudes(isou,ifreq);
            end
        
            if ifreq == 1
                timingdata(2) = etime(clock,t00);
                t00 = clock;
            end
            
            % End of CASE 1, PART 2 
            % Q_firstterm has been computed
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % CASE 1,PART 3: Starting from Q_firstterm, iterate to find Qfinal.
            %
            % Two for-loops: iterations (iiter) and submatrices (nn)

            if showtext >= 2
                disp('      Starting to iterate')
            end
        
            Q = Q_firstterm;
            Qfinal = Q;   
            
            for iiter = 1:(controlparameters.difforder-2)
                if showtext >= 3
                   disp(['         Iteration no. ',int2str(iiter),' of ',int2str(controlparameters.difforder-2)])
                end
                Q_ackum = Q*0;              
                shortlistprev = 0;
                Hsubdata = [];
                ivuse = [];

                for nn = 1:size(Hsubmatrixdata.listofsubmatrices,1)
                    ivoutvec = Hsubmatrixdata.bigmatrixstartnums(Hsubmatrixdata.listofsubmatrices(nn,2)):Hsubmatrixdata.bigmatrixendnums(Hsubmatrixdata.listofsubmatrices(nn,2));
                    ivinvec  = Hsubmatrixdata.bigmatrixstartnums(Hsubmatrixdata.listofsubmatrices(nn,3)):Hsubmatrixdata.bigmatrixendnums(Hsubmatrixdata.listofsubmatrices(nn,3));
                      
                    shortlistnumber = Hsubmatrixdata.reftoshortlist(nn);
                    if shortlistnumber ~= shortlistprev
                           Hsubdata = Hsubdatalists{nn}; 
                           ivuse = ivuselists{nn}; 
                           nHsub1 = length(ivoutvec);
                           row_Hsub = rem(ivuse-1, nHsub1) + 1;
                           col_Hsub = (ivuse-row_Hsub)/nHsub1 + 1;
                           shortlistprev = shortlistnumber;                   
                    end
                    
                    Q_ackum(1:ivoutvec(end)) = Q_ackum(1:ivoutvec(end)) + accumarray(ivoutvec(1)-1+row_Hsub,Hsubdata.*Q(ivinvec(1)-1+col_Hsub));                      
                end
                Q = Q_ackum; 
                Qfinal = Qfinal + Q; 
            end
                    
            doesQsegmenthavevalues = zeros(length(Hsubmatrixdata.edgepairlist(:,1)),1);
            for iicounter = 1:length(Hsubmatrixdata.edgepairlist(:,1))
               sectionofbigmatrix = (Hsubmatrixdata.bigmatrixstartnums(iicounter):Hsubmatrixdata.bigmatrixendnums(iicounter)); 
               doesQsegmenthavevalues(iicounter) = any(Qfinal(sectionofbigmatrix));
            end        

            if ifreq == 1
                timingdata(3) = etime(clock,t00);
                t00 = clock;
            end

            % End of CASE 1, PART 3 
            % Qfinal has been computed
            % Possibly save to a file
            
            if filehandlingparameters.saveinteqsousigs == 1
               filename_sousigs = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_f',int2str(ifreq),'_sousigs.mat'];  
               eval(['save ''',filename_sousigs,''' Qfinal Q_firstterm doesQsegmenthavevalues'])             
            end

        else    % This 'else' part implies that controlparameters.loadinteqsousig = 1
            filename_sousigs = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_f',int2str(ifreq),'_sousigs.mat'];  
            eval(['load ''',filename_sousigs,'''']) 
			disp(['      Using existing edgesourcesignalsfile: ',filename_sousigs])
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CASE 1, PART 4: Qfinal has been computed, or loaded from a file
        %
        % Propagate Qfinal to the external receiver(s)
        
        for irec = 1:nreceivers
            if showtext >=2
                disp(['      Receiver number ',int2str(irec)]) 
            end
            
            [Fmatrix,ivproblematic] = EDcalcpropagatematrix(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,isthinhole,Rdata.vispartedgesfromr(:,irec),...
                    frequency,controlparameters.Rstart,...
                    Rdata.rRsho(Rdata.reftoshortlistR(:,irec)),Rdata.thetaRsho(Rdata.reftoshortlistR(:,irec)),Rdata.zRsho(Rdata.reftoshortlistR(:,irec)),doesQsegmenthavevalues,filehandlingparameters.showtext);
                
             P_receiver(ifreq,irec) = Fmatrix*Qfinal;
             if filehandlingparameters.savediff2result == 1
                extraoutputdata.tfinteqdiff_nodiff2(ifreq,irec) = Fmatrix*(Qfinal-Q_firstterm);
             end
             
            if ifreq == 1
                timingdata(4) = etime(clock,t00);
            end
        end % for irec = 1:nreceivers
    
        % End of CASE 1, PART 4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CASE 2: doaddsources = 0  and we have several sources
    % We might also have the case that doallSRcombinations = 0.
    % 
    % PART 2: Create Q_firstterm for each source
    % There will be an outer loop over sources
         
    elseif Sdata.doaddsources == 0   % This means we have several sources, and they should not be added
        if ifreq==1
            P_receiver = zeros(nfrequencies,nreceivers,nsources);
        end
        for isou = 1:nsources
            
            if filehandlingparameters.loadinteqsousigs == 0
                
                ISOU = int2str(isou); 
    
                if showtext >= 2
                    disp(['Source no. ',ISOU,' of ',int2str(nsources)])
                end
                
                if strcmp(Sdata.sourcetype,'monopole') == 1
                    Q_firstterm = EDinteg_souterm(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,...
                        controlparameters,Sdata.vispartedgesfroms(:,isou),Sdata.vispartedgesfroms_start(:,isou),...
                        Sdata.vispartedgesfroms_end(:,isou),frequency,gaussvectors,Sdata.rSsho(Sdata.reftoshortlistS(:,isou)),...
                        Sdata.thetaSsho(Sdata.reftoshortlistS(:,isou)),Sdata.zSsho(Sdata.reftoshortlistS(:,isou)),filehandlingparameters.showtext);
                elseif strcmp(Sdata.sourcetype,'polygonpiston') == 1
                    Q_firstterm_sumpistongaussorder = 0;
                    for np = 1:size(Sdata.pistongausscoordinates{1},1)
                            Q_firstterm_addition = EDinteg_souterm(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,...                    
                            controlparameters,Sdata.vispartedgesfroms(:,isou),Sdata.vispartedgesfroms_start(:,isou),...
                            Sdata.vispartedgesfroms_end(:,isou),frequency,gaussvectors,Sdata.rSsho(Sdata.reftoshortlistS(:,isou,np)),...
                            Sdata.thetaSsho(Sdata.reftoshortlistS(:,isou,np)),Sdata.zSsho(Sdata.reftoshortlistS(:,isou,np)),filehandlingparameters.showtext);
                        Q_firstterm_sumpistongaussorder = Q_firstterm_sumpistongaussorder + Q_firstterm_addition*squeeze(Sdata.pistongaussweights{isou}(np,1));
                    end
                    Q_firstterm = Q_firstterm_sumpistongaussorder;

                end

                Q_firstterm = Q_firstterm*Sdata.sourceamplitudes(isou,ifreq);
                              
                % End of CASE 2, PART 2 
                % Q_firstterm has been computed
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % CASE 2,PART 3: Starting from Q_firstterm, iterate to find Qfinal.
                %
                % Two for-loops: iterations (iiter) and submatrices (nn)
                  
                if showtext >= 2
                    disp('      Starting to iterate')
                end
                Q = Q_firstterm;
                Qfinal = Q;   
                
                for iiter = 1:(controlparameters.difforder-2)
                   if showtext >= 3
                       disp(['         Iteration no. ',int2str(iiter),' of ',int2str(controlparameters.difforder-2)])
                   end
                   Q_ackum = Q*0;              
                   shortlistprev = 0;
                   Hsubdata = [];
                   ivuse = [];
                   for nn = 1:size(Hsubmatrixdata.listofsubmatrices,1)
                        ivoutvec = Hsubmatrixdata.bigmatrixstartnums(Hsubmatrixdata.listofsubmatrices(nn,2)):Hsubmatrixdata.bigmatrixendnums(Hsubmatrixdata.listofsubmatrices(nn,2));
                        ivinvec  = Hsubmatrixdata.bigmatrixstartnums(Hsubmatrixdata.listofsubmatrices(nn,3)):Hsubmatrixdata.bigmatrixendnums(Hsubmatrixdata.listofsubmatrices(nn,3));

                        shortlistnumber = Hsubmatrixdata.reftoshortlist(nn);
                        if shortlistnumber ~= shortlistprev
                           Hsubdata = Hsubdatalists{nn}; 
                           ivuse = ivuselists{nn}; 
                           nHsub1 = length(ivoutvec);
                           row_Hsub = rem(ivuse-1, nHsub1) + 1;
                           col_Hsub = (ivuse-row_Hsub)/nHsub1 + 1;
                           shortlistprev = shortlistnumber;                   
                        end
                         Q_ackum(1:ivoutvec(end)) = Q_ackum(1:ivoutvec(end)) + accumarray(ivoutvec(1)-1+row_Hsub,Hsubdata.*Q(ivinvec(1)-1+col_Hsub));                      
                    end
                    Q = Q_ackum; 
                    Qfinal = Qfinal + Q; 
                end % for iiter = ....
                       
                doesQsegmenthavevalues = zeros(length(Hsubmatrixdata.edgepairlist(:,1)),1);
                for iicounter = 1:length(Hsubmatrixdata.edgepairlist(:,1))
                   sectionofbigmatrix = (Hsubmatrixdata.bigmatrixstartnums(iicounter):Hsubmatrixdata.bigmatrixendnums(iicounter)); 
                   doesQsegmenthavevalues(iicounter) = any(Qfinal(sectionofbigmatrix));
                end

                % End of CASE 2, PART 3 
                % Qfinal has been computed
                % Possibly save to a file

                if filehandlingparameters.saveinteqsousigs == 1
                   filename_sousigs = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_',ISOU,'_f',int2str(ifreq),'_sousigs.mat'];  
                   eval(['save ''',filename_sousigs,''' Qfinal Q_firstterm doesQsegmenthavevalues'])             
                end
            
            else    % This 'else' part implies that filehandlingparameters.loadinteqsousig = 1                
                filename_sousigs = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_',ISOU,'_f',int2str(ifreq),'_sousigs.mat'];  
                eval(['load ''',filename_sousigs,''''])              	
                disp(['   Using existing edgesourcesignalsfile: ',filename_sousigs])
            end  

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % CASE 2, PART 4: Qfinal has been computed, or loaded from a file
            %
            % Propagate Qfinal to the external receiver(s)
            
%             if inteqsettings.handlepropagation < 3
            
            if Sdata.doallSRcombinations == 1
                for irec = 1:nreceivers

                     if showtext >=2
                        disp(['Receiver number ',int2str(irec)]) 
                     end
                     [Fmatrix,ivproblematic] = EDcalcpropagatematrix(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,isthinhole,...
                        Rdata.vispartedgesfromr(:,irec),frequency,controlparameters.Rstart,...
                        Rdata.rRsho(Rdata.reftoshortlistR(:,irec)),...
                        Rdata.thetaRsho(Rdata.reftoshortlistR(:,irec)),Rdata.zRsho(Rdata.reftoshortlistR(:,irec)),doesQsegmenthavevalues,filehandlingparameters.showtext);

                    if ifreq == 1 && ~isempty(ivproblematic)
                        existinglength = size(ivproblematicmatrix,2);
                        newlength = length(ivproblematic);
                        if newlength > existinglength
                            ivproblematicmatrix = [ivproblematicmatrix zeros(nreceivers,newlength-existinglength)];
                        end
                        ivproblematicmatrix(irec,1:newlength) = ivproblematic;
                    end

                     P_receiver(ifreq,irec,isou) = Fmatrix*Qfinal;
                     if filehandlingparameters.savediff2result == 1
                         extraoutputdata.tfinteqdiff_nodiff2(ifreq,irec,isou) = Fmatrix*(Qfinal-Q_firstterm);
                     end

                end %  for irec = 1:nreceivers
            else
                irec = isou;
                  if showtext >=2
                    disp(['Receiver number ',int2str(irec)]) 
                 end
                 [Fmatrix,ivproblematic] = EDcalcpropagatematrix(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,isthinhole,...
                    Rdata.vispartedgesfromr(:,irec),frequency,controlparameters.Rstart,...
                    Rdata.rRsho(Rdata.reftoshortlistR(:,irec)),...
                    Rdata.thetaRsho(Rdata.reftoshortlistR(:,irec)),Rdata.zRsho(Rdata.reftoshortlistR(:,irec)),doesQsegmenthavevalues,filehandlingparameters.showtext);

                if ifreq == 1 && ~isempty(ivproblematic)
                    existinglength = size(ivproblematicmatrix,2);
                    newlength = length(ivproblematic);
                    if newlength > existinglength
                        ivproblematicmatrix = [ivproblematicmatrix zeros(nreceivers,newlength-existinglength)];
                    end
                    ivproblematicmatrix(irec,1:newlength) = ivproblematic;
                end

                 P_receiver(ifreq,irec,isou) = Fmatrix*Qfinal;
                 if filehandlingparameters.savediff2result == 1
                     extraoutputdata.tfinteqdiff_nodiff2(ifreq,irec,isou) = Fmatrix*(Qfinal-Q_firstterm);
                 end
               
            end
        end   % for isou = 1:nsources

    end  % if Sdata.doaddsources == 1 || nsources == 1
end  % for ifreq = 1:nfrequencies     

elapsedtimehodtf = etime(clock,t00);
tfinteqdiff = P_receiver;

EDinputdatahash = EDinputdatahash_tfinteq;

eval(['save(''',desiredname,''',''tfinteqdiff'',''extraoutputdata'',''EDinputdatahash'',''timingdata'',''elapsedtimehodtf'');'])

    
