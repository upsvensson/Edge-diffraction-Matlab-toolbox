function [geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters] = ...
    EDcheckinputstructs(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters)
% EDcheckinputstructs checks the input data structs to EDmain_convex
% and sets default values.
% 
% Input parameters:
%   geoinputdata     .geoinputfile        (specified, or a file open
%                                             window is given - unless the fields .corners and .planecorners are given)
%                   .corners             (optional; alternative to
%                                             .geoinputfile. But, the use of this option
%                                             requires that filehandlingparameters.outputdirectory and .filestem are specified)
%                   .planecorners        (optional; alternative to
%                                             .geoinputfile)
%                   .planerefltypes      (optional; gives the reflection
%                                         factors of the planes: 1,0,-1.
%                                         default: ones(nplanes,1) )
%                   .firstcornertoskip   (optional, default: 1e6)
%                   .listofcornerstoskip (optional, default [])
%                   .freefieldcase       (optional, default: 0)
%                   .planeseesplanestrategy  (optional, defaul: 0)
%   Sinputdata         .coordinates         (obligatory)
%                   .doaddsources        (default: 0 = no)
%                   .sourceamplitudes    Only used if doaddsources = 1
%                                        (default:
%                                         ones(nsources,nfrequencies))
%                   .doallSRcombinations  (default: 1 = yes)
%                   .sourcetype           'monopole' (default) or
%                                         'polygonpiston'
%                   .pistoncornercoordinates (obligatory if .sourcetype =
%                   'polygonpiston'). Matrix size[npistoncorners,3].
%                   .pistoncornernumbers   (obligatory if .sourcetype =
%                   'polygonpiston'). Matrix size [npistons,maxncornersperpiston].
%                   .pistonplanes          (obligatory if .sourcetype =
%                   'polygonpiston'). List of size [npistons,1].
%                   .nedgesubs           (optional, default 2)
%   Rinputdata         .coordinates         (obligatory)
%                   .nedgesubs           (optional, default 2)
%   envdata         .cair                (default: 344)
%                   .rhoair              (default: 1.21)
%   controlparameters   .fs              (default: 44100) Only relevant if
%                                        controlparameters.docalcir = 1.
%                   .directsound         (default: 1 = yes)
%                   .difforder           (default: 10)
%                   .savealldifforders   (default: 0) Only relevant if
%                                        controlparameters.docalcir = 1.
%                   .saveindividualfirstdiff (default: 0) Only relevant if
%                                        controlparameters.docalcir = 1.
%                   .skipfirstorder      (default: 0)
%             Maximum one of .docalctf,.docalcir,docalctf_ESIEBEM can be 1.
%                   .docalctf            (default: 0) Only on 
%                   .docalcir            (default: 0)
%                   .docalctf_ESIEBEM    (default: 0)
%                   .Rstart              (default: 0)
%                   .frequencies         (obligatory if .docalctf = 1 or
%                                        docalctf_ESIEBEM = 1)
%                   .discretizationtype  (default: 2 = G-L)
%                   .ngauss              (default: 16, but ignored if
%                                        .docalcir = 1)
%                   .surfacegaussorder   (default: 5, but ignored if 
%                                         .docalctf_ESIEBEM = 0)
%                   .HODirelemsize       (default:
%                                       2*2.^(-[0:controlparameters.difforder-1]))
%   filehandlingparameters    .outputdirectory  (default: the folder of the geoinputfile)  
%                   .filestem        (default: name of the cad-file, with an underscore + running integer)
%                   .savecadgeofile      (default: 0)
%                   .saveSRdatafiles     (default: 1)
%                   .saveeddatafile      (default: 1)
%                   .saveed2datafile      (default: 1)
%                   .savesubmatrixdata    (default: 1)
%                   .saveinteqsousigs     (default: 0)
%                   .loadinteqsousigs     (default: 0)
%                   .savepathsfile        (default: 1)
%                   .saveISEStree         (default: 0) Not used by
%                                         EDmain_convexESIE
%                   .savelogfile          (default: 1)
%                   .savediff2result      (default: 0)
%                   .savehodpaths         (default: 0) Used only if
%                                         .docalcir = 1
%                   .showtext             (default: 1)
% 
% Peter Svensson 22 March 2024 (peter.svensson@ntnu.no)
% 
% [geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters] = ...
% EDcheckinputstructs(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters);

% 24 Nov. 2017 First version
% 28 Nov. 2017 Cleaned code a bit
% 29 Nov. 2017 Corrected mistake: saveSRdatafiles was called
%              saveSRinputdatafiles. Adjusted to the new indata option:
%              specified corners and planecorners matrices instead of a CAD
%              file.
% 30 Nov. 2017 Fixed a bug where docalctf, instead of docalcir, was set to
%              0. Changed defaults to create eddata and SRfiles, because of
%              plotting the model. Changed one field from logfilename to
%              savelogfile. Removed the field lineending.
% 13 Dec. 2017 Added the input field sourceamplitudes
% 12 Jan. 2018 Forced doaddsources to be 1, if the number of sources = 1.
%              Also followed one yellow recommendation: numel instead of
%              prod.
% 17 Jan 2018  Added the input field planecornertype
% 18 Jan 2018  Fixed a bug with the sourceamplitudes; they were forced to
% one by mistake. Changed default for savelogfile to 1.
% 22 Jan 2018 Removed the field planecornertype in the struct geoinputdata
% 22 Jan 2018 Moved some controlparameter fields away from the convex TF
% case: fs, nedgepoints_visibility, docalcir, specorder.
% 22 Jan 2018 Changed the defaults for saving files.
% 26 Jan 2018 Changed the defaults for saving files. Also changed the
% default output directory to include "results" in the path.
% 26 Jan 2018 Introduced Sinputdata.doallSRcombinations. Default value 1.
% 28 Jan 2018 First version for EDmain_convexESIE_ir
% 31 Jan 2018 Corrected the handling of sourceamplitudes (the freq.
% dependence was not implemented correctly.
% 6 Feb 2018 Introduced a check if number of receivers/sources was zero
% 8 Feb 2018 v 0.109 Introduced a new parameter:
% controlparameters.skipfirstorder (default = 0).
% 14 Feb 2018 v0.112 Small change, assigning a default value to filestem if
% it was not given.
% 14 Feb 2018 Added a check if the functions 'DataHash.m' and 'lgwt.m' are
% available.
% 14 Feb 2018 Added a test if the source and receiver coordinates has the
% right number of columns (3).
% 15 Feb 2018 Removed the savesetupfile. Introduced the saveed2datafile,
% default = 1. Changed defaults for savepathsfile and savesubmatrixdata to
% 1. Stopped assigning default value to saveISEStree. Introduced the
% parameter suppressresultrecycling, default = 0.
% 1 Mar 2018 Introduced the version 3: EDmain_convexESIEBEM
% 1 Mar 2018 Stripped away irrelevant fields; otherwise the DataHash
% doesn't always recognize what is identical.
% 15 Mar 2018 Added the field filehandlingparameters.savehodpaths
% 16 Mar 2018 Combined two error messages
% 21 Mar 2018 Made a few changes to the controlparameters: gave default
% values for docalctf and docalcir. Removed specorder. Introduced the
% .savealldifforders parameter. 
% 21 Mar 2018 Changed savealldifforders to controlparameters instead of
% filehandlingparameters
% 7 Apr 2018 Introduced controlparameters.saveindividualfirstdiff
% 21 Apr 2018 Removed the extra "results" directory in the output
% directory.
% 28 May 2018 Added the input field geoinputdata.planerefltypes
% 22 May 2019 Fixed an error with sourceamplitudes; they did not have the
% right orientation.
% 3 June 2020 Fixed a but with sourceamplitudes that was found by EDdebug:
% If the user had specified a constant sourceamplitudes, it wasn't expanded
% to a matrix of size [nsources, nfrequencies].
% 14 March 2021 The section "% Check the struct Rinputdata" was moved to
% before "% Check the struct Sinputdata"
% 28 Sep. 2023 Adapted to EDmain_convex which does both tf and ir
% 2 Oct. 2023 Changed so that only one of docalftf, docalcir, docalctf_ESIEBEM
% can be set to zero.
% 4 Oct. 2023 Changed default difforder to 10
% 12 Oct. 2023 Introduced the field geoinputdata.freefieldcase
% 26 Oct. 2023 Introduced new fields in Sinputdata which define piston
% sources.
% 29 Oct. 2023 Introduced the new fields .listofcornerstoskip and 
% .planeseesplanestrategy in geoinputdata. Introduced nedgesubs for Sdata
% and Rdata. Added the new field .HODirelemsize
% 30 Oct. 2023 Introduced the piston source checking
% 4 Nov. 2023 Calculates the piston areas.
% 22 Nov. 2023 Fixed a bug: regular polygons had not got the gaussweights
% as a vector, only as a single value.
% 29 Nov. 2023 Changed the name of the filed pistongausspoints to
% pistongaussorder.
% 22 March 2024 Updated the function header (there were references to
% EDmain_convexESIE etc, which was removed).
% 22 March 2024 Removed the requirement that one of the calculation options
% must be 1. If no calculation option is set to 1, then only geometry files
% will be generated, and that is very useful when the geometry is being
% developed.
% 22 March 2024 Added the field listofedgestoskip to geoinputdata

% if nargin < 7
%     disp('ERROR: the input parameter EDmaincase was not specified')
%     return
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if the needed non-EDtoolbox functions are available.

file1ok = (exist('DataHash.m','file') == 2);
file2ok = (exist('lgwt.m','file') == 2);

if file1ok == 0 && file2ok == 0
   error('ERROR: Matlab can not find the functions lgwt.m and DataHash.m. Please download them from Mathworks') 
end
if file1ok == 0
   error('ERROR: Matlab can not find the function DataHash.m. Please download it from Mathworks') 
end
if file2ok == 0
   error('ERROR: Matlab can not find the function lgwt.m. Please download it from Mathworks') 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct geoinputdata

if ~isstruct(geoinputdata)

	[CADfile,CADfilepath] = uigetfile('*.*','Please select the cadfile');
    [~,CADfile,~] = fileparts([CADfilepath,CADfile]);

    CADfile = [CADfilepath,CADfile];
    geoinputdata = struct ('geoinputfile',CADfile);

    [infilepath,CADfilestem] = fileparts(geoinputdata.geoinputfile);
    if ~isfield(filehandlingparameters,'outputdirectory')
        filehandlingparameters.outputdirectory = [infilepath,filesep,'results'];
    end   
    if ~isfield(filehandlingparameters,'filestem')
        filehandlingparameters.filestem = CADfilestem;
    end
end
if ~isfield(geoinputdata,'freefieldcase')
    geoinputdata.freefieldcase = 0;
end
if geoinputdata.freefieldcase == 1
    geoinputdata.corners = [];
    geoinputdata.planecorners = [];
else
    if isfield(geoinputdata,'geoinputfile')
        [infilepath,CADfilestem] = fileparts(geoinputdata.geoinputfile);
        if ~isfield(filehandlingparameters,'outputdirectory')
            filehandlingparameters.outputdirectory = [infilepath,filesep,'results'];
        end    
        if ~isfield(filehandlingparameters,'filestem')
            filehandlingparameters.filestem = CADfilestem;
        end
    else
    
        if ~isfield(geoinputdata,'corners') || ~isfield(geoinputdata,'planecorners')
    	    [CADfile,CADfilepath] = uigetfile('*.*','Please select the cadfile');
            [~,CADfile,~] = fileparts([CADfilepath,CADfile]);
    
            CADfile = [CADfilepath,CADfile];
            geoinputdata.geoinputfile = CADfile;
            [infilepath,CADfilestem] = fileparts(geoinputdata.geoinputfile);        
            if ~isfield(filehandlingparameters,'outputdirectory')
                filehandlingparameters.outputdirectory = [infilepath,filesep,'results'];
            end    
        else
            if isfield(filehandlingparameters,'outputdirectory') == 0 || isfield(filehandlingparameters,'filestem') == 0
                error('ERROR: When you give the geometry input in the form of data matrices, you must specify filehandlingparameters.outputdirectory and .filestem')            
            end
            if ~isfield(geoinputdata,'planerefltypes')
               nplanes = size(geoinputdata.planecorners,1);
               geoinputdata.planerefltypes = ones(nplanes,1);
            end
        end
    end
    if ~isfield(geoinputdata,'firstcornertoskip')
        geoinputdata.firstcornertoskip = 1e6;
    end
    if ~isfield(geoinputdata,'listofcornerstoskip')
        geoinputdata.listofcornerstoskip = [];
    end
    if ~isfield(geoinputdata,'listofedgestoskip')
        geoinputdata.listofedgestoskip = [];
    else
        ncols = size(geoinputdata.listofedgestoskip,2);
        if ncols > 2
            error('ERROR: geoinputdata.listofedgestoskip has more than two columns')
        else
            if any(any(geoinputdata.listofedgestoskip==0))
                error('ERROR: geoinputdata.listofedgestoskip contains some zeros')
            end
        end
    end
    if ~isfield(geoinputdata,'planeseesplanestrategy')
        geoinputdata.planeseesplanestrategy = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct Rinputdata

if ~isstruct(Rinputdata)
    error('ERROR 1: receiver coordinates were not specified')
end
if ~isfield(Rinputdata,'coordinates')
    error('ERROR 2: receiver coordinates were not specified')
end
nreceivers = size(Rinputdata.coordinates,1);
if nreceivers == 0
     error('ERROR 3: receiver coordinates were not specified')            
end
ncolumns = size(Rinputdata.coordinates,2);
if ncolumns ~= 3
   error(['ERROR: check your receiver coordinates; there were ',int2str(ncolumns),' columns rather than 3']) 
end
if ~isfield(Rinputdata,'nedgesubs')
    Rinputdata.nedgesubs = 2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct Sinputdata, but the field .sourceamplitudes is handled
% further down, since the values depend on whether it is a TD case or an FD
% case.

if ~isstruct(Sinputdata)
    error('ERROR: the source was not specified')
end
% sourcetype is not specified, the default is monopole
if ~isfield(Sinputdata,'sourcetype')
    Sinputdata.sourcetype = 'monopole';
    Sinputdata.pistoncornercoordinates = [];
    Sinputdata.pistoncornernumbers = [];
    Sinputdata.pistonplanes = [];    
end
if strcmp(Sinputdata.sourcetype,'monopole')
    if ~isfield(Sinputdata,'coordinates') 
        error('ERROR: you must specify monopole source coordinates')
    end
    [nsources,ncolumns] = size(Sinputdata.coordinates);
    if nsources == 0
         error('ERROR: source coordinates were not specified')            
    end
    if ncolumns ~= 3
       error(['ERROR: check your sources coordinates; there were ',int2str(ncolumns),' columns rather than 3']) 
    end
    Sinputdata.pistoncornercoordinates = [];
    Sinputdata.pistoncornernumbers = [];
    Sinputdata.pistonplanes = []; 
    Sinputdata.pistongaussorder = 0;
    
elseif strcmp(Sinputdata.sourcetype,'polygonpiston')
   if ~isfield(Sinputdata,'pistoncornercoordinates') || ~isfield(Sinputdata,'pistoncornernumbers') || ~isfield(Sinputdata,'pistonplanes')
        error('ERROR: when you specify a polygonpiston sources, you must specify corner coordinates and numbers, and piston plane(s).')
    end
    npistons = size(Sinputdata.pistoncornernumbers,1);
    [npoints,ncolumns] = size(Sinputdata.pistoncornercoordinates);
    if npistons == 0
        error('ERROR: You must define at least one piston source by specifying the corner numbers')
    else
        nsources = npistons;
    end
    if npoints < 3
        error('ERROR: You must define at least three pistoncornercoordinates')
    end
    if ncolumns ~= 3
        error('ERROR: There should be 3 columns for the pistoncornercoordinates')
    end        
    if ~isfield(Sinputdata,'pistongaussorder')
        Sinputdata.pistongaussorder = 3;
    end

    pistonmidpointcoordinates = zeros(npistons,3);   
    ncornersperpiston = sum(Sinputdata.pistoncornernumbers.'>0).';
    pistonareavec = zeros(npistons,1);

    Sinputdata.pistongausscoordinates = cell(npistons,1);
    Sinputdata.pistongaussweights = cell(npistons,1);

    for ii = 1:npistons
        corners_onepiston = Sinputdata.pistoncornercoordinates...
            (Sinputdata.pistoncornernumbers(ii,1:ncornersperpiston(ii)),:);
        pistonmidpointcoordinates(ii,:) = mean(corners_onepiston);

        % Calculate the area of each piston. Implemented for triangles
        % and rectangles, using a version of Heron's formula for the
        % triangle: A = 1/4*sqrt( 4*a^2*b^2 - (a^2 + b^2 - c^2)^2 )
        % where a,b,c = triangle side lengths
        corners_onepiston_expanded = [corners_onepiston;corners_onepiston(1,:)];
        sidevectors = diff(corners_onepiston_expanded);
        sidelengthsquared = sum(sidevectors.^2,2);

        if ncornersperpiston(ii) == 3
            pistonarea = 1/4*sqrt( 4*sidelengthsquared(1)*sidelengthsquared(2) ...
                - (sidelengthsquared(1) + sidelengthsquared(2) - sidelengthsquared(3))^2  );
            switch Sinputdata.pistongaussorder
                case {0,1}
                   gaussvalues = [1/3 1/3 1]; 
                case {2}
                    gaussvalues = [2/3 1/6 1/3;...
                                   1/6 2/3 1/3;
                                   1/6 1/6 1/3];
                case {3}
                    gaussvalues = [1/3 1/3 -9/32; ...
                        0.2 0.6 25/96; ...
                        0.2 0.2 25/96; ...
                        0.6 0.2 25/96];
                case {4}
                    gaussvalues = [0.108103018168070  0.445948490915965 0.223381589678011; ...
                                   0.445948490915965  0.108103018168070 0.223381589678011; ...
                                   0.445948490915965  0.445948490915965 0.223381589678011; ...
                                   0.816847572980459  0.091576213509771 0.109951743655322; ...
                                   0.091576213509771  0.816847572980459 0.109951743655322; ...
                                   0.091576213509771  0.091576213509771 0.109951743655322];
                case {5}   % 7 gauss points
                    gaussvalues = [1/3  1/3  0.225;...
                                   0.059715871789770 0.470142064105115 0.132394152788506;...
                                   0.470142064105115 0.059715871789770 0.132394152788506;...
                                   0.470142064105115 0.470142064105115 0.132394152788506;...
                                   0.797426985353087 0.101286507323456 0.125939180544827;...
                                   0.101286507323456 0.797426985353087 0.125939180544827;...
                                   0.101286507323456 0.101286507323456 0.125939180544827];
                case{6}  % 12 gauss points
                    gaussvalues = [0.501426509658179 0.249286745170910 0.116786275726379;...
                                   0.249286745170910 0.501426509658179 0.116786275726379;...
                                   0.249286745170910 0.249286745170910 0.116786275726379;...
                                   0.873821971016996 0.063089014491502 0.050844906370207;...
                                   0.063089014491502 0.873821971016996 0.050844906370207;...
                                   0.063089014491502 0.063089014491502 0.050844906370207;...                           
                                   0.053145049844817 0.310352451033784 0.082851075618374;...
                                   0.310352451033784 0.053145049844817 0.082851075618374;...
                                   0.053145049844817 0.636502499121399 0.082851075618374;...
                                   0.636502499121399 0.053145049844817 0.082851075618374;...
                                   0.636502499121399 0.310352451033784 0.082851075618374;...
                                   0.310352451033784 0.636502499121399 0.082851075618374];
                case{7,8} % 16 gauss points 
                      gaussvalues = [1/3 1/3 0.144315607677787;...
                                     0.081414823414554 0.459292588292723 0.095091634267285;...
                                     0.459292588292723 0.081414823414554  0.095091634267285;...
                                     0.459292588292723 0.459292588292723 0.095091634267285;...
                                     0.658861384496480 0.170569307751760 0.103217370534718;...
                                     0.170569307751760 0.658861384496480 0.103217370534718;...
                                     0.170569307751760 0.170569307751760 0.103217370534718;...
                                     0.898905543365938 0.050547228317031 0.032458497623198;...
                                     0.050547228317031 0.898905543365938 0.032458497623198;...
                                     0.050547228317031 0.050547228317031 0.032458497623198;...
                                     0.008394777409958 0.263112829634638 0.027230314174435;...
                                     0.728492392955404 0.008394777409958 0.027230314174435;...
                                     0.008394777409958 0.728492392955404 0.027230314174435;...
                                     0.728492392955404 0.263112829634638 0.027230314174435;...
                                     0.263112829634638 0.008394777409958 0.027230314174435;...
                                     0.263112829634638 0.728492392955404 0.027230314174435];
                                                               
                otherwise
                    error(['ERROR: triangle quadrature not implemented for this number: ',int2str(surfacegaussorder)])
            end
            n2 = size(gaussvalues,1);  
%             corners_x = planedata.corners(planedata.planecorners(planenumber,1:3),1);
%             corners_y = planedata.corners(planedata.planecorners(planenumber,1:3),2);
%             corners_z = planedata.corners(planedata.planecorners(planenumber,1:3),3);
            corners_x = corners_onepiston(:,1);
            corners_y = corners_onepiston(:,2);
            corners_z = corners_onepiston(:,3);
            
            % Convert the normalized triangle coordinates to the 3D coordinates
            
            psi=[gaussvalues(:,1) gaussvalues(:,2) 1-gaussvalues(:,1)-gaussvalues(:,2)];
            
            xgausspoints = psi*corners_x;
            ygausspoints = psi*corners_y;
            zgausspoints = psi*corners_z;
            Sinputdata.pistongausscoordinates{ii} = [xgausspoints ygausspoints zgausspoints];

            Sinputdata.pistongaussweights{ii} = gaussvalues(:,3);

        elseif ncornersperpiston(ii) == 4
            extradistsquared = sum( (corners_onepiston_expanded(2,:) - corners_onepiston_expanded(4,:)).^2 );
            triarea1 =  1/4*sqrt( 4*sidelengthsquared(1)*sidelengthsquared(2) ...
                - (sidelengthsquared(1) + sidelengthsquared(2) -extradistsquared)^2  );
            triarea2 =  1/4*sqrt( 4*sidelengthsquared(3)*sidelengthsquared(4) ...
                - (sidelengthsquared(3) + sidelengthsquared(3) -extradistsquared)^2  );
            pistonarea = triarea1 + triarea2;

           n2 = Sinputdata.pistongaussorder^2;
           [x,w] = lgwt(Sinputdata.pistongaussorder,0,1);
           
           x = x(end:-1:1);
           xy = [repmat(x,Sinputdata.pistongaussorder,1) reshape(x(:,ones(1,Sinputdata.pistongaussorder)).',n2,1)];
           ww = prod([repmat(w,Sinputdata.pistongaussorder,1) reshape(w(:,ones(1,Sinputdata.pistongaussorder)).',n2,1)],2);   
            
           c1 = corners_onepiston(1,:);
           c2 = corners_onepiston(2,:);
           c3 = corners_onepiston(3,:);
           c4 = corners_onepiston(4,:);
           xvec1 = c2 - c1;
           alen = norm(xvec1);
           xvec2 = c3 - c4;
           clen = norm(xvec2);
           yvec1 = c4 - c1;
           dlen = norm(yvec1);
           yvec2 = c3 - c2;
           blen = norm(yvec2);
           crossvec1 = c3 - c1;
           plen = norm(crossvec1);
           crossvec2 = c4 - c2;
           qlen = norm(crossvec2);         

           startpointsx = c1(ones(n2,1),:);
           startpointsx = startpointsx + xvec1(ones(n2,1),:).*xy(:,[1 1 1]);
    
           endpointsx = c4(ones(n2,1),:);
           endpointsx = endpointsx + xvec2(ones(n2,1),:).*xy(:,[1 1 1]);
           
           startpointsy = c1(ones(n2,1),:);
           startpointsy = startpointsy + yvec1(ones(n2,1),:).*xy(:,[2 2 2]);
    
           endpointsy = c2(ones(n2,1),:);
           endpointsy = endpointsy + yvec2(ones(n2,1),:).*xy(:,[2 2 2]);

           startpointsx = c1(ones(n2,1),:);
           startpointsx = startpointsx + xvec1(ones(n2,1),:).*xy(:,[1 1 1]);
    
           endpointsx = c4(ones(n2,1),:);
           endpointsx = endpointsx + xvec2(ones(n2,1),:).*xy(:,[1 1 1]);
           
           startpointsy = c1(ones(n2,1),:);
           startpointsy = startpointsy + yvec1(ones(n2,1),:).*xy(:,[2 2 2]);
    
           endpointsy = c2(ones(n2,1),:);
           endpointsy = endpointsy + yvec2(ones(n2,1),:).*xy(:,[2 2 2]);
    
           Sinputdata.pistongausscoordinates{ii} = startpointsx + (endpointsx - startpointsx).*xy(:,[2 2 2]);
           Sinputdata.pistongaussweights{ii} = ww;
        else
            % Check if the polygon is regular; then we can compute its area
            if std(sidelengthsquared)/mean(sidelengthsquared) > 1e-9
                error('ERROR: piston areas can only be calculated for 3- and 4-sided pistons, and regular polygons.')
            else
                cornerradius = mean( EDcalcdist(corners_onepiston,pistonmidpointcoordinates(ii,:)) );
                pistonarea = ncornersperpiston(ii)/2*cornerradius^2*sin(2*pi/ncornersperpiston(ii));

                % We should not use cornerradius; it would be a bit more
                % accurate to use the theoretical circle radius
                circlecoordinates = GEOcirclepoints(cornerradius,Sinputdata.pistongaussorder);
                
                Sinputdata.pistongausscoordinates{ii} = circlecoordinates + pistonmidpointcoordinates(ii,:)
                Sinputdata.pistongaussweights{ii} = 1/size(circlecoordinates,1)*ones(size(circlecoordinates,1),1);
                disp('WARNING! Regular polygon piston has only been properly implemented in the plane with nvec = [0 0 1]')
            end
        end
        pistonareavec(ii) = pistonarea;
    end
    Sinputdata.coordinates = pistonmidpointcoordinates;
    Sinputdata.ncornersperpiston = ncornersperpiston;
    Sinputdata.pistonareas = pistonareavec;

    if ~isfield(Sinputdata,'pistonplanes')
        error('ERROR: When Sinputdata.sourcetype is polygonpiston, then the field .pistonplanes must be specified.')
    else
        if size(Sinputdata.pistonplanes(:),1) ~= npistons
            error('ERROR: the fields .pistoncornernumbers and .pistonplanes must have the same number of rows.')
        end
    end

else
    error('ERROR: Only monopole and polygonpistonsources have been defined')
end

if ~isfield(Sinputdata,'doaddsources')
    Sinputdata.doaddsources = 0;
end
if nsources == 1
    Sinputdata.doaddsources = 1;
end    
if ~isfield(Sinputdata,'doallSRcombinations')
    Sinputdata.doallSRcombinations = 1;
else
    if Sinputdata.doallSRcombinations == 0 && nsources~=nreceivers
        disp(['   nsources = ',int2str(nsources)])
        disp(['   nreceivers = ',int2str(nreceivers)])
        
       error('ERROR: doallSRcombinations was set to 0, but the number of sources was not the same as the number of receivers'); 
    end
end

if ~isfield(Sinputdata,'nedgesubs')
    Sinputdata.nedgesubs = 2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct envdata

if ~isstruct(envdata)
    error('ERROR: the struct envdata was not specified')
end
if ~isfield(envdata,'cair')
    envdata.cair = 344;
end
if ~isfield(envdata,'rhoair')
    envdata.rhoair = 1.21;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct controlparameters

if ~isstruct(controlparameters)
     error('ERROR: the struct controlparameters was not specified')           
end

% First some general parameters 

if ~isfield(controlparameters,'directsound')
    controlparameters.directsound = 1;
end
if ~isfield(controlparameters,'skipfirstorder')
    controlparameters.skipfirstorder = 0;
end
if ~isfield(controlparameters,'Rstart')
    controlparameters.Rstart = 0;
end
if ~isfield(controlparameters,'difforder') 
    if geoinputdata.freefieldcase == 1
        controlparameters.difforder = 0;
    else
        disp('   INFO: controlparameters.difforder was not specified. It is given the value 10.')
        controlparameters.difforder = 10;
    end
end

% Then other parameters 

if ~isfield(controlparameters,'docalctf')
    controlparameters.docalctf = 0;
%    disp('WARNING: controlparameters.docalctf was not specified. It is given the value 0.')
end
if ~isfield(controlparameters,'docalcir')
    controlparameters.docalcir = 0;
%    disp('WARNING: controlparameters.docalcir was not specified. It is given the value 0.')
end
if (controlparameters.docalctf == 1) && (controlparameters.docalcir == 1)
    error('EROR: docalcir and docalctf can not both be set to 1')
end        
if ~isfield(controlparameters,'docalctf_ESIEBEM')
    controlparameters.docalctf_ESIEBEM = 0;
%    disp('WARNING: controlparameters.docalctf_ESIEBEM was not specified. It is given the value 0.')
else
    if controlparameters.docalctf_ESIEBEM == 1 
        if geoinputdata.freefieldcase == 1
            error('ERROR: controlparameters.docalctf_ESIEBEM and geoinputdata.freefieldcase can not both be set to 1')
        end
        if controlparameters.docalctf == 1
            error('ERROR: controlparameters.docalctf and controlparameters.docalctf_ESIEBEM can not both be set to 1')
        end
        if controlparameters.docalcir == 1
            error('ERROR: controlparameters.docalcir and controlparameters.docalctf_ESIEBEM can not both be set to 1')
        end
    end
end
% if (controlparameters.docalctf == 0) && (controlparameters.docalcir == 0) && (controlparameters.docalctf_ESIEBEM == 0)
%     error('ERROR: You must choose a calulcation method: one of (docalcir, docalctf, docalctf_ESIEBEM) must be 1')
% end

if controlparameters.docalctf == 1 || controlparameters.docalctf_ESIEBEM == 1
    if ~isfield(controlparameters,'frequencies')
        error('ERROR: controlparameters.frequencies were not specified')
    end
    if ~isfield(controlparameters,'ngauss')
        disp('   INFO: controlparameters.ngauss wasnt set; it is given the default value 16.')
        controlparameters.ngauss = 16;
    end
    if ~isfield(controlparameters,'discretizationtype')
%        disp('   INFO: controlparameters.discretizationtype wasnt set; it is given the default value 2.')
        controlparameters.discretizationtype = 2;
    end
else
    controlparameters.frequencies = [];
    controlparameters.ngauss = 0;
    controlparameters.discretizationtype = 0;
end
nfrequencies = length(controlparameters.frequencies);

if controlparameters.docalctf_ESIEBEM == 1
    if ~isfield(controlparameters,'surfacegaussorder')
        disp('   INFO: controlparameters.surfacegaussorder wasnt set; it is given the default value 5.')
        controlparameters.surfacegaussorder = 5;
    end
else
    controlparameters.surfacegaussorder = 0;
end

if controlparameters.docalcir == 1
    if ~isfield(controlparameters,'fs')
        disp('   INFO: controlparameters.fs wasn''t set; it is given the default value 44100.')
        controlparameters.fs = 44100;
    end   
    if ~isfield(controlparameters,'HODelemsize')
        controlparameters.HODelemsize = 2*2.^(-[0:controlparameters.difforder-1]);
    end   
else
   controlparameters.fs = 44100;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the field .sourceamplitudes of the struct Sinputdata

sourceamplitudes_default = 0;

% We need to treat the frequency-domain cases differently from the
% time-domain cases. But, if .sourceamplitudes was not specified, it gets
% the default value 1 for all sources (and frequencies).

if ~isfield(Sinputdata,'sourceamplitudes')
%     Sinputdata.sourceamplitudes = ones(nsources,1);    
    Sinputdata.sourceamplitudes = 1;    
    sourceamplitudes_default = 1;
end

[n1,n2] = size(Sinputdata.sourceamplitudes);

if n1 > 1
    if n1 ~= nsources
        error(['ERROR: The Sinputdata.sourceamplitudes input parameter must have one row, or one row per source, but it had the size ',int2str(n1),' by ',int2str(n2)])
    else % So, we know that there is one amplitude per source. Check if there is one value per freq.
        if controlparameters.docalctf == 1 || controlparameters.docalctf_ESIEBEM == 1      
            if n2 == 1
                Sinputdata.sourceamplitudes = Sinputdata.sourceamplitudes(:,ones(1,nfrequencies));
            else
                if n2 ~= nfrequencies
                    error(['ERROR: The Sinputdata.sourceamplitudes input parameter must have one column, or one column per frequency, but it had the size ',int2str(n1),' by ',int2str(n2)])
                end
            end
        else
            if n2 > 1
               error(['ERROR: The Sinputdata.sourceamplitudes input parameter must, for TD calculations, have one column, but it had the size ',int2str(n1),' by ',int2str(n2)])
            end
        end
    end
else % Here we know that there is just one row. Check if there is one value per freq. (or a single value)
    if n2 > 1
        if controlparameters.docalctf == 1 || controlparameters.docalctf_ESIEBEM == 1      
            if n2 ~= nfrequencies
                error(['ERROR: The Sinputdata.sourceamplitudes input parameter must have one column, or one column per frequency, but it had the size ',int2str(n1),' by ',int2str(n2)])
            else % Here we know that n2 == nfrequencies. Expand to number of sources
                Sinputdata.sourceamplitudes = Sinputdata.sourceamplitudes(ones(nsources,1),:);
            end
        else
            error(['ERROR: The Sinputdata.sourceamplitudes input parameter must, for TD calculations, have one column, but it had the size ',int2str(n1),' by ',int2str(n2)])            
        end
    else % Here we know that n1 = 1 and n2 = 1
        if controlparameters.docalctf == 1 || controlparameters.docalctf_ESIEBEM == 1      
            Sinputdata.sourceamplitudes = Sinputdata.sourceamplitudes(ones(nsources,1),ones(1,nfrequencies));
        else
            Sinputdata.sourceamplitudes = Sinputdata.sourceamplitudes(ones(nsources,1));            
        end
    end
end

if controlparameters.docalcir == 1
    if ~isfield(controlparameters,'savealldifforders')
        controlparameters.savealldifforders = 0;
    end
    if ~isfield(controlparameters,'saveindividualfirstdiff')
        controlparameters.saveindividualfirstdiff = 0;
    end
else
    controlparameters.savealldifforders = 0;
    controlparameters.saveindividualfirstdiff = 0;
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct filehandlingparameters
    
if exist([filehandlingparameters.outputdirectory],'dir') ~=7
      mkdir([filehandlingparameters.outputdirectory])
end

if ~isfield(filehandlingparameters,'filestem')
    if exist('CADfilestem','var') == 1
        ifile = 1;
        filename1 = [filehandlingparameters.outputdirectory,filesep,'results',filesep,CADfilestem,'_',int2str(ifile),'_tf.mat'];
        filename2 = [filehandlingparameters.outputdirectory,filesep,'results',filesep,CADfilestem,'_',int2str(ifile),'_ir.mat'];
        if exist(filename1,'file') == 2 || exist(filename2,'file') == 2
            failedtofindfile = 0;
            while (failedtofindfile == 0)
                ifile = ifile + 1;
                filename1 = [filehandlingparameters.outputdirectory,filesep,'results',filesep,CADfilestem,'_',int2str(ifile),'_tf.mat'];
                filename2 = [filehandlingparameters.outputdirectory,filesep,'results',filesep,CADfilestem,'_',int2str(ifile),'_ir.mat'];
                if exist(filename1,'file') ~= 2 && exist(filename2,'file') ~= 2
                    failedtofindfile = 1;
                end
            end
        end
        filehandlingparameters.filestem = [CADfilestem,'_',int2str(ifile)];
    else
       filehandlingparameters.filestem = Filestem;        
    end
end
% if ~isfield(filehandlingparameters,'savesetupfile')
%     filehandlingparameters.savesetupfile = 1;
% end
if ~isfield(filehandlingparameters,'suppressresultrecycling')
    filehandlingparameters.suppressresultrecycling = 0;
end
if ~isfield(filehandlingparameters,'showtext')
    filehandlingparameters.showtext = 1;
end
if ~isfield(filehandlingparameters,'savecadgeofile')
    filehandlingparameters.savecadgeofile = 0;
end
if ~isfield(filehandlingparameters,'saveSRdatafiles')
    filehandlingparameters.saveSRdatafiles = 1;
end
if ~isfield(filehandlingparameters,'saveeddatafile')
    filehandlingparameters.saveeddatafile = 1;
end
if ~isfield(filehandlingparameters,'saveed2datafile')
    filehandlingparameters.saveed2datafile = 1;
end
if ~isfield(filehandlingparameters,'savesubmatrixdata')
    filehandlingparameters.savesubmatrixdata = 1;
end
if ~isfield(filehandlingparameters,'saveinteqsousigs')
    filehandlingparameters.saveinteqsousigs = 0;
end
if ~isfield(filehandlingparameters,'loadinteqsousigs')
    filehandlingparameters.loadinteqsousigs = 0;
end
if ~isfield(filehandlingparameters,'savepathsfile')
    filehandlingparameters.savepathsfile = 1;
end

if controlparameters.docalcir == 1
    if ~isfield(filehandlingparameters,'saveISEStree')
        filehandlingparameters.saveISEStree = 0;
    end
    if ~isfield(filehandlingparameters,'savehodpaths')
        filehandlingparameters.savehodpaths = 0;
    end
else
    filehandlingparameters.savehodpaths = 0;
    filehandlingparameters.saveISEStree = 0;
end
if ~isfield(filehandlingparameters,'savediff2result')
    filehandlingparameters.savediff2result = 0;
end
if ~isfield(filehandlingparameters,'savelogfile')
    filehandlingparameters.savelogfile = 1;
end

