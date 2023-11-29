function dirsoundok = EDdirectsound(corners,planecorners,planeeqs,planenvecs,...
    S,R,visplanesfroms,visplanesfromr,canplaneobstruct,minvals,maxvals,ncornersperplanevec,showtext)
% EDdirectsound - Checks if the direct sound is valid.
% Checks if the direct sound path from a number of source
% points to a number of receiver points are obstructed by any planes.
%
% Input parameters:
%   eddatafile      The name of the file that contains all edge related data.
%                   This file will be loaded.
%   S, R, visplanesfroms, visplanesfromr
%                   Data that should have been passed directly from the
%                   srdatafile.
%
% Output parameter:
%   dirsoundok      1 or 0, telling if the sound path is unobstructed or not.
%
% Peter Svensson (peter.svensson@ntnu.no) 16 March 2021
%
% Uses the function EDchkISvisible
%
% dirsoundok = EDdirectsound(corners,planecorners,planeeqs,planenvecs,...
% S,R,visplanesfroms,visplanesfromr,canplaneobstruct,minvals,maxvals,ncornersperplanevec);

% 29 Sept. 2003 Functioning version
% 28 Nov. 2017 Copied to EDtoolbox. Introduced the non-global showtext
% input parameter
% 16 Mar 2021 Adapted to change in EDchkISvisible.

% global showtext
% showtext = 1;

% ncorners = size(corners,1);
% nplanes = size(planecorners,1);

%----------------------------------------------------------
% 	Check if the direct sound is obscured
%
%   Pick out the planes that are possible. Only potentially obstructing planes are necessary to check. 
%   Planes that are seen by the source or by the receiver, but not both,
%   are potentially obstructing.

%planesareseen = full( double(visplanesfroms)==1 + double(visplanesfromr)==1  );
planesareseen = full( (double(visplanesfroms)==0 | double(visplanesfromr)==0) & (double(visplanesfroms)+double(visplanesfromr)~=0)  );

planestocheck = find(planesareseen.*double(canplaneobstruct).');
nplanestocheck = length(planestocheck);
onesvec = ones(nplanestocheck,1);

% [hitplanes,reflpoints,edgehits,edgehitpoints,cornerhits,cornerhitpoints] = EDchkISvisible(S(onesvec,:),R(onesvec,:),planeeqs(planestocheck,4),planenvecs(planestocheck,:),minvals(planestocheck,:),...
%    maxvals(planestocheck,:),planecorners(planestocheck,:),corners,ncornersperplanevec(planestocheck));
[hitplanes,~,~,~,~,~,~,~] = EDchkISvisible(S(onesvec,:),R(onesvec,:),planeeqs(planestocheck,4),planenvecs(planestocheck,:),minvals(planestocheck,:),...
   maxvals(planestocheck,:),planecorners(planestocheck,:),corners,ncornersperplanevec(planestocheck));

if isempty(hitplanes)
	dirsoundok = 1;
	if showtext >= 3
		disp('         direct sound OK')
	end
else
	dirsoundok = 0;
    if showtext >= 3
        disp('         direct sound obscured')
        if showtext >= 4
            printvec = int2str(planestocheck(hitplanes(1)));
            for iiprint = 2:length(hitplanes)
   		        printvec = [printvec,' ',int2str(planestocheck(hitplanes(iiprint)))];      
            end
            disp('          by:')
  		    disp(['         ',printvec])
        end
    end
end
