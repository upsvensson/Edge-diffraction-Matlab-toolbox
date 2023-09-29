function [Hsubmatrixdata,outpar,existingfilename] = EDinteg_submatrixstructure(edgelengthvec,closwedangvec,...
    inteq_ngauss,inteq_discretizationtype,edgetoedgedata,planesatedge,EDversionnumber,inpar)
% EDinteg_submatrixstructure - Develop the edge-to-edge matrix structure for the ED int.eq.
% 
% A version 1 of this function was used up to v 0.221 of EDtoolbox
% and version 2 after that.
%
% Input parameters:
% 	edgelengthvec       A list [nedge,1] of the lenghts of each edge.
% 	closwedangvec       A list [nedge,1] of the close-wedge angle of each
%                       edge.
%   inteq_ngauss                No. of disc. points along the longest edge
%   inteq_discretizationtype    [0,1,2]
%   edgetoedgedata      struct
%   planesatedge        A matrix, size [nedges,nplanes] 
%   EDversionnumber
%   inpar               This parameter is different for version 1 and
%                       version 2 of this function. 
%       v1 inpar should be showtext (optional)
%                       0 -> no text displayed. Value > 0 -> text displayed.
%       v2 inpar should be filehandlingparameters (obligatory)
%                       filehandlingparameters is a struct which
%                       contains the field showtext.
%
% Output parameters:
%   Hsubmatrixdata      Struct with the fields
%       Hsubmatrix
%       listofsubmatrices
%       edgepairlist
%       edgetripletlist
%       submatrixcounter
%       nuniquesubmatrices
%       nedgeelems
%       bigmatrixstartnums
%       bigmatrixendnums
%       shortlist
%       reftoshortlist
%       isthinplanetriplet
%   outpar              This parameter is different for version 1 and
%                       version 2 of this function. 
%       v1 outpar = EDinputdatahash       
%                       This is a string of characters which is
%                       uniquely representing the input data.
%                       An existing result file with the same value of this
%                       EDinputdatahash will be reused.
%       v2 outpar = elapsedtimesubmatrix
%                       This tells how long time was used inside this
%                       function. If an existing file was reused, then
%                       elapsedtimesubmatrix has a second value which tells
%                       how much time was used for the existing file.
%   existingfilename    For v2 of this function, if an existing file was
%                       found that was reused, then the reused file
%                       name is given here. If no existing file could be 
%                       reused then this variable is empty. For v1 of
%                       this function, this variable is also empty.%
%
% Uses the functions  EDdistelements, EDrecycleresultfiles from EDtoolbox
% Uses the function DataHash from Matlab Central
%
% Peter Svensson (peter.svensson@ntnu.no) 28 Sep. 2023 
%
% [Hsubmatrixdata,outpar,existingfilename] = EDinteg_submatrixstructure...
%    (edgelengthvec,closwedangvec,inteq_ngauss,inteq_discretizationtype,...
%    edgetoedgedata,planesatedge,EDversionnumber,inpar)

% 29-1-2013 Cleaned up quite much and avoided looping through all the
%           empty combinations. Exported much more data. Introduced
%           edgepairlist and edgetripletlist
% 30-1-2013 Introduced the listofsubmatrices and isthinplanetriplet output parameters.
%           Introduced a reshuffling of the rows in edgepairlist, and in
%           the Hsubmatrices such that sizes change less often.
% 6-2-2014  Removed a temporary change, so now, the
%           number of edge elements depends on the length of each edge.
% 2 Nov. 2017 Changed the display of text somewhat
% 27 Nov. 2017 Copied from ESIE2toolbox. Moved the edge discretization part
%               out of this function.
%  28 Nov. 2017 Introduced the non-global showtext input parameter
% 29 Nov. 2017 Changed call from ESIE2distelements to EDdistelements
% 5 Dec. 2017 Allow odd numbers of edge points. Stored one extra field:
%               nuniquesubmatrices
% 15 Dec. 2017 Turned off the code which suppresses symmetrycompressions
% 25 Jan 2018 Made sure there are at least 2 edge elements per edge.
% 8 Feb 2018 Introduced the EDinputdatahash
% 28 Sep. 2023 Implemented version 2 of this function while maintaining
% compatibility with the old "version 1". v2 moves the check if an existing
% file can be reused inside this function. Also updated load and save to
% the function call form, which avoids problems with spaces in file names.

t00 = clock;

if nargin < 8  % Must be the old version
	functionversion = 1;
	showtext = 0;
else % nargin = 8 -> could be the old or new version
	if isstruct(inpar)
		functionversion = 2;
		filehandlingparameters = inpar;
        showtext = filehandlingparameters.showtext;
	else
		functionversion = 1;
		showtext = inpar;
	end
end

EDinputdatastruct = struct('edgelengthvec',edgelengthvec,'closwedangvec',closwedangvec,...
    'inteq_ngauss',inteq_ngauss,'inteq_discretizationtype',inteq_discretizationtype,...
    'edgetoedgedata',edgetoedgedata,'planesatedge',planesatedge,'EDversionnumber',EDversionnumber);
EDinputdatahash = DataHash(EDinputdatastruct);

if functionversion == 1
	outpar = EDinputdatahash;
	existingfilename = '';
elseif functionversion == 2

    %---------------------------------------------------------------
    % Sort out the file business: can an existing file be used?
    % Then copy the existing file to a new copy. Should the data be saved in a file? 
    
    if filehandlingparameters.suppressresultrecycling == 1
        foundmatch = 0;
        existingfilename = '';
    else
        [foundmatch,existingfilename] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_submatrixdata',EDinputdatahash);
    end
    
    desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_submatrixdata.mat'];
    
    if foundmatch == 1
        eval(['load(''',existingfilename,''')'])
        if ~strcmp(existingfilename,desiredname)
            copyfile(existingfilename,desiredname);
        end
        elapsedtimesubmatrix_new = etime(clock,t00);
        elapsedtimesubmatrix = [elapsedtimesubmatrix_new elapsedtimesubmatrix];
        outpar = elapsedtimesubmatrix;
        return
    end
end

% ----------------------------------------------------------------------------------------------
% Create the two-index list for the edge-source column matrix.
% Check the edgeseespartialedge if the edge-to-edge combination should be
% included or not.

symmetrycompression = 1;
if showtext >= 3
    disp(['      EDinteg_e2estructure: symmetrycompression = ',int2str(symmetrycompression)])
end

nedges = size(edgetoedgedata.edgeseespartialedge,1);

iv = (1:nedges).';

% iv1 and iv2 will be the edge numbers
%   iv1 (first column) is "to-edge"
%   iv2 (second column) is "from-edge"

iv1 = repmat(iv,nedges,1);
iv2 = reshape(repmat(iv.',nedges,1),nedges*nedges,1);

iv = [iv1 iv2];
ivremove = iv(:,1)==iv(:,2);
iv(ivremove,:) = [];

% Remove the combinations where an edge can not see another edge
for ii = 1:size(iv,1)
    if edgetoedgedata.edgeseespartialedge(iv(ii,1),iv(ii,2)) == 0
        iv(ii,:) = [0 0];
    end
end
ivkeep = (iv(:,1)~=0);
iv = iv(ivkeep,:);

nedgepairs = size(iv,1);

if nedges < 256
   edgepairlist = uint8(iv);
   edgetripletlist = zeros(nedgepairs^2,3,'uint8');
elseif nedges < 65536
   edgepairlist = uint16(iv);
   edgetripletlist = zeros(nedgepairs^2,3,'uint16');
end

% activeedges = unique([edgepairlist(:,1);edgepairlist(:,2)]);

% ----------------------------------------------------------------------------------------------
% Change the order of the rows in edgepairlist such that the changes in
% vector size happen less often.
%
% First the number of discretization points is determined for each edge.

% Before 5 Dec. 2017: only even numbers of edge points were used
%  nedgeelems = ceil(inteq_ngauss*(edgelengthvec/max(edgelengthvec))/2)*2;
nedgeelems = ceil(inteq_ngauss*(edgelengthvec/max(edgelengthvec)));

% 25 Jan 2018: Make sure there are at least 2 edge elements per edge.
iv = (nedgeelems == 1);
nedgeelems(iv) = 2;

[~,sortvec] = sortrows([nedgeelems(edgepairlist(:,1)) nedgeelems(edgepairlist(:,2))]);
edgepairlist = [edgepairlist(sortvec,1) edgepairlist(sortvec,2)];

% ----------------------------------------------------------------------------------------------
% Make two lists, bigmatrixstartnums and bigmatrixendnums, which for each
% edge 2 -> edge 1 combo specifies which rows in the Q-vector are used.

bigmatrixelemvec = zeros(nedgepairs,1);
for ii = 1:nedgepairs
    bigmatrixelemvec(ii) = nedgeelems(edgepairlist(ii,1))*nedgeelems(edgepairlist(ii,2));
end
% nbig = sum(bigmatrixelemvec);
bigmatrixendnums = cumsum(bigmatrixelemvec);
bigmatrixstartnums = [1;bigmatrixendnums(1:end-1)+1];

% ----------------------------------------------------------------------------------------------
% Prepare for the big H-matrix by first creating a matrix, Hsubmatrix,
% showing which submatrices should be included.
%
% This is determined by the fact that the combination
%   Source   edge A -> edge B   
%   Target   edge C -> edge D
% is only valid if edge B == edge C and then the submatrix corresponds to:
%   From edge A, via edge B/C, to edge D
%
% Give those submatrices a unique integer, counting up from 1.
%
% edgetripletlist will be a matrix which for row i contains the edge-triplet
%   to-edge, via-edge, from-edge
% for submatrix number i

Hsubmatrix = zeros(nedgepairs,nedgepairs,'uint32');

submatrixcounter = uint32(0);

for ii = 1:nedgepairs
    ivec = find( edgepairlist(ii,2) == edgepairlist(:,1) );
    nivec = uint32(length(ivec));
    Hsubmatrix(ii,ivec) = submatrixcounter + (1:nivec);
    
    edgetripletlist(submatrixcounter+(1:nivec),:) = [ double(edgepairlist(ii,1))*ones(nivec,1) double(edgepairlist(ii,2))*ones(nivec,1)  edgepairlist(ivec,2)];
    submatrixcounter = submatrixcounter + nivec;
end

edgetripletlist = edgetripletlist(1:submatrixcounter,:);

% ----------------------------------------------------------------------------------------------
% Change the order of the submatrices such that they change as little as
% possible from row to row. This is done either 
% 1. using "symmetrycompression", or 
% 2. Based on number of edge elements.

% 1. Compute the edge-to-cylindrical coordinates and find if there are
% identical combos (within 1e-5 accuracy). All edge-related coordinate
% values are used.
%
% The "reftoshortlist" will be a matrix which for each submatrix gives the
% number in a "shortlist" which will restore the correct values. The
% "shortlist" then contains the uniques sets of edge-related coordinate
% values.
% 
% 2. Sort based on the number of edge elements for edge 2 and 3 in
% edgetriplets.

if symmetrycompression == 1
    if showtext >= 3
        disp('      Looking for identical edge-to-edge-to-edge combos')
    end

    refto = edgetoedgedata.reftoshortlistE;

    datamatrix = zeros(submatrixcounter,12);
    alledge3 = edgetripletlist(:,1);
    alledge2 = edgetripletlist(:,2);
    alledge1 = edgetripletlist(:,3);
    nyvec = pi./(2*pi-closwedangvec);
    allny = nyvec(alledge2);
    alldz = edgelengthvec(alledge2)./nedgeelems(alledge2); 
    datamatrix(:,1:3) = [allny alldz edgelengthvec(alledge2)];
    
    for ii = 1:submatrixcounter        
        datamatrix(ii,4) = abs( edgetoedgedata.thetae1sho(refto(alledge1(ii),alledge2(ii))) - edgetoedgedata.thetae1sho(refto(alledge3(ii),alledge2(ii))) );
        datamatrix(ii,5:12) = [edgetoedgedata.ze1sho(refto(alledge3(ii),alledge2(ii))),edgetoedgedata.ze2sho(refto(alledge3(ii),alledge2(ii))),...
                                   edgetoedgedata.re1sho(refto(alledge3(ii),alledge2(ii))),edgetoedgedata.re2sho(refto(alledge3(ii),alledge2(ii))),...
                                   edgetoedgedata.ze1sho(refto(alledge1(ii),alledge2(ii))),edgetoedgedata.ze2sho(refto(alledge1(ii),alledge2(ii))),...
                                   edgetoedgedata.re1sho(refto(alledge1(ii),alledge2(ii))),edgetoedgedata.re2sho(refto(alledge1(ii),alledge2(ii)))];
    end
    
    [shortlist,~,reftoshortlist] = unique(round(datamatrix*1e5)/1e5,'rows');
    nuniquesubmatrices = size(shortlist,1);    

%     if size(shortlist,1) > round(length(reftoshortlist)*0.8)
%         symmetrycompression = -1;
%     end
        
    if showtext >= 3
        disp(' ')
        if symmetrycompression == 1
            disp(['      Can use symmetry in Hsub-matrices: will compute and save ',int2str(size(shortlist,1)),' instead of ',int2str(length(reftoshortlist)),' sub-matrices'])
        else
            disp(['      Looked for symmetry in Hsub-matrices but could not find any. Will compute and save ',int2str(length(reftoshortlist)),' sub-matrices'])            
        end
    end
else
    reftoshortlist = [];
end

listofsubmatrices = find(Hsubmatrix);
[~,listsortvec] = sort(Hsubmatrix(listofsubmatrices));
listofsubmatrices = listofsubmatrices(listsortvec);

if symmetrycompression < 1
    [~,sortvec] = sortrows(nedgeelems(edgetripletlist(:,2:3)));
    Hsubmatrix(listofsubmatrices(sortvec)) = (1:submatrixcounter);
    edgetripletlist = edgetripletlist(sortvec,:);
    if ~isempty(reftoshortlist)
        reftoshortlist = reftoshortlist(sortvec);
    end
    listofsubmatrices = listofsubmatrices(sortvec);
    reftoshortlist = [];
else
   [~,reftosortvec] = sort(reftoshortlist);
   Hsubmatrix(listofsubmatrices(reftosortvec)) = (1:submatrixcounter);
    edgetripletlist = edgetripletlist(reftosortvec,:);
    reftoshortlist = reftoshortlist(reftosortvec);
    listofsubmatrices = listofsubmatrices(reftosortvec);
end

jjlist = floor((listofsubmatrices-1)/nedgepairs)+1;
iilist = listofsubmatrices-(floor((listofsubmatrices-1)/nedgepairs))*nedgepairs;

listofsubmatrices = [listofsubmatrices iilist jjlist];

% ----------------------------------------------------------------------------------------------
% Figure out which edge triplets and which edge pairs correspond to thin planes

% First some special cases

thinplanecheckready = 0;

iv = find(closwedangvec == 0);

if isempty(iv)
    % No thin edges in entire model
    isthinplanetriplet  = zeros(submatrixcounter,1,'int8');
    isthinplaneedgepair = zeros(nedgepairs,1,'int8');
    thinplanecheckready = 1;
else
    if length(iv) == length(closwedangvec)
        % Only thin edges - then we have to check if it is a disc.
        allplanes = reshape(planesatedge,size(planesatedge,1)*size(planesatedge,2),1);
        if length(unique(allplanes)) == 2
            isthinplanetriplet  = ones(submatrixcounter,1);
            isthinplaneedgepair = ones(nedgepairs,1);
            thinplanecheckready = 1;
        end
    end
end

if thinplanecheckready == 0
    involvedplanes = [planesatedge(edgepairlist(:,1),:) planesatedge(edgepairlist(:,2),:)];
    diffinvolved = diff(sort(involvedplanes.'));
    isthinplaneedgepair = int8(sum(diffinvolved>0)==1);
    
    involvedplanes = [planesatedge(edgetripletlist(:,1),:) planesatedge(edgetripletlist(:,2),:) planesatedge(edgetripletlist(:,3),:)];
    diffinvolved = diff(sort(involvedplanes.'));
    isthinplanetriplet = int8(sum(diffinvolved>0)==1);
end

% ----------------------------------------------------------------------------------------------
% Run through the Gauss-Legendre discretization for each number of edge
% points that occur, and store the data in a matrix.

uniquelistofedgepoints = unique(nedgeelems);

quadraturematrix_pos     = sparse(zeros(max(uniquelistofedgepoints),max(uniquelistofedgepoints)));
quadraturematrix_weights = sparse(zeros(max(uniquelistofedgepoints),max(uniquelistofedgepoints)));

for ii = 1:length(uniquelistofedgepoints)
    n = uniquelistofedgepoints(ii);
    [n1vec,weightvec] = EDdistelements(n,inteq_discretizationtype);
    quadraturematrix_pos(n,1:n)     = n1vec(:).';
    quadraturematrix_weights(n,1:n) = weightvec(:).';
end

Hsubmatrixdata = struct('Hsubmatrix',Hsubmatrix);
Hsubmatrixdata.listofsubmatrices = listofsubmatrices;
Hsubmatrixdata.edgepairlist = edgepairlist;
Hsubmatrixdata.edgetripletlist = edgetripletlist;
Hsubmatrixdata.submatrixcounter = submatrixcounter;
Hsubmatrixdata.nuniquesubmatrices = nuniquesubmatrices;
Hsubmatrixdata.nedgeelems = nedgeelems;
Hsubmatrixdata.bigmatrixstartnums = bigmatrixstartnums;
Hsubmatrixdata.bigmatrixendnums = bigmatrixendnums;
Hsubmatrixdata.reftoshortlist = reftoshortlist;
Hsubmatrixdata.isthinplanetriplet = isthinplanetriplet;
Hsubmatrixdata.isthinplaneedgepair = isthinplaneedgepair;
Hsubmatrixdata.quadraturematrix_pos = quadraturematrix_pos;
Hsubmatrixdata.quadraturematrix_weights = quadraturematrix_weights;

if functionversion == 2
    elapsedtimesubmatrix = etime(clock,t00);
    outpar = elapsedtimesubmatrix;
    
    if filehandlingparameters.savesubmatrixdata == 1
    	eval(['save(''',desiredname,''',''Hsubmatrixdata'',''EDinputdatahash'',''elapsedtimesubmatrix'');'])
    end
end








