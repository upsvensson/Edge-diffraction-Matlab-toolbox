% EDdebug3 runs EDmain_convex_time for several defined tests, to check its handling of
% different numbers of sources, receivers, and frequencies
% as well as values of difforder = 0,1,2, doaddsources = 0,1
% sourceamplitudes = 1 or ones(nfrequencies,nsources)
% 
% Peter Svensson 13 Apr 2018 (peter.svensson@ntnu.no)

% 31 Jan 2018 First version
% 7 Mar 2018 Version for the ESIEBEM
% 11 Apr 2018 Version for the EDmain_convex_time
% 13 Apr 2018 Introduced the saveindividualfirstdiff as varied input
% parameter

[EDversionnumber,changedate,changetime] = EDgetversion;

% Two sources, three receivers, four frequencies
% doaddsources = 0 or 1

% sources = [0.04 0 0.00001;0 0 -0.32001];
% sources = [0.04 0 0.1;0 0 0.1;-0.04 0 0.1];%0 0 -0.32001];
sources = [1 2 3;4 -5 -6;-10 -12 13];
receivers = [0 0 -2; 1 1 1];
% frequencies = [1 10 20 40];

sources = EDcirclepoints(0.1,2);
sources(:,3) = 0.00001;

% sources = sources(1:2,:);

mfile = mfilename('fullpath');
[infilepath,filestem] = fileparts(mfile);

% corners = [     -0.2000   -0.4400   -0.3200
%     0.2000   -0.4400   -0.3200
%     0.2000    0.2000   -0.3200
%    -0.2000    0.2000   -0.3200
%    -0.2000   -0.4400         0
%     0.2000   -0.4400         0
%     0.2000    0.2000         0
%    -0.2000    0.2000         0];
corners = [     -0.2000   -0.4   -0.32
    0.2000   -0.4   -0.32
    0.2000    0.2000   -0.32
   -0.2000    0.2000   -0.32
   -0.2000   -0.4         0
    0.2000   -0.4         0
    0.2000    0.2000         0
   -0.2000    0.2000         0];

planecorners = [   1     4     3     2
     5     6     7     8
     1     2     6     5
     3     4     8     7
     2     3     7     6
     1     5     8     4];

geoinputdata = struct('corners',corners,'planecorners',planecorners);

filehandlingparameters = struct('outputdirectory',[infilepath,filesep,'results']);
filehandlingparameters.filestem = filestem;
filehandlingparameters.savelogfile = 0;
filehandlingparameters.showtext = 0;
controlparameters = struct('ngauss',8);
controlparameters.discretizationtype = 2;
controlparameters.Rstart = 0;
controlparameters.fs = 8000;
% controlparameters.surfacegaussorder = 2;
envdata.cair = 344;
Sinputdata = struct;

souvar = [1 size(sources,1)];
recvar = [1 size(receivers,1)];
addsouvar = [0 1];
diffordervar = [0 1 2];
souamp = [1 2];

casecounter = 0;
for ii = 1:length(souvar)
    Sinputdata.coordinates = sources(1:souvar(ii),:);
    nsources = size(Sinputdata.coordinates,1);
    for jj = 1:length(recvar)
        Rinputdata.coordinates = receivers(1:recvar(jj),:);
        for kk = [ 1 2]
            controlparameters.saveindividualfirstdiff = kk-1;
            for ll = 1:length(diffordervar)
                controlparameters.difforder = diffordervar(ll);
                for mm = 1:length(addsouvar)
                    Sinputdata.doaddsources = addsouvar(mm);   
                    for nn = 1:length(souamp)
                        if nn == 1
                            Sinputdata.sourceamplitudes = 1;
                        else
                           Sinputdata.sourceamplitudes = ones(1,nsources); 
                        end
                        
                        casecounter = casecounter + 1;
                        disp(['Case no. ',int2str(casecounter),': ',int2str(ii),' ',int2str(jj),' ',int2str(ll),' ',int2str(mm),' ',int2str(nn)])
                        EDmain_convex_time(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);
                        eval(['load ',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_ir.mat irdirect irdiff'])

                        if controlparameters.saveindividualfirstdiff == 0
                            if any(any(any(irdiff))) && controlparameters.difforder == 0
                                error('ERROR: difforder was set to zero but irdiff got some result')
                            end
                            if ~any(any(any(irdiff))) && controlparameters.difforder > 0
                                disp('WARNING: difforder was set > 0 but irdiff got no result. Check if this is as expected')
                            end
                        else
                            if controlparameters.difforder == 0
                                if any(any(any(irdiff))) && controlparameters.difforder == 0
                                    error('ERROR: difforder was set to zero but irdiff got some result')
                                end                                
                            else
                                if ~any(any(any(irdiff{1,1}.irvectors)))
                                    disp('WARNING: difforder was set > 0 but irdiff got no result. Check if this is as expected')
                                end  
                            end
                            
                                                      
                        end
                        if controlparameters.difforder > 1
                            eval(['load ',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_irhod.mat irhod'])
                            if any(any(any(irhod{1,1}))) && controlparameters.difforder <2
                                error('ERROR: difforder was set <2 but irhod got some result')
                            end
                            if ~any(any(any(irhod{1,1}))) && controlparameters.difforder > 1
                                error('ERROR: difforder was set > 1 but irhod got no result')
                            end
                        end
                    end
                end                

            end
        end
    end
end


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Test 1: results should be 4 * 3 * 2
% % % 
% % % No diffraction at all
% % 
% % 
% %     Sinputdata.doaddsources = 0;
% %     controlparameters.difforder = 0;
% %     
% %     EDmain_convexESIE(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);
% % 
% %     eval(['load ',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tf.mat tfdirect tfdiff'])
% % 
% %     if ndims(tfdirect) < 3
% %         error('Test 1: ERROR 1')
% %     else
% %         if size(tfdirect,1) ~= 4
% %            error('Test 1: ERROR 2') 
% %         end
% %         if size(tfdirect,2) ~= 3
% %            error('Test 1: ERROR 3') 
% %         end
% %     end
% %     
% %     
% %     
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Test 2: results should be 4 * 3 
% % % 
% % % difforder = 1
% % 
% % 
% %     Sinputdata = struct('coordinates',sources);
% %     Sinputdata.doaddsources = 1;
% %     Rinputdata = struct('coordinates',receivers);
% %     controlparameters.frequencies = frequencies;
% %     controlparameters.difforder = 1;
% %     
% %     EDmain_convexESIE(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);
% % 
% %     eval(['load ',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tf.mat tfdirect tfdiff'])
% % 
% %     if ndims(tfdirect) > 2
% %         error('Test 2: ERROR 1')
% %     else
% %         if size(tfdirect,1) ~= 4
% %            error('Test 2: ERROR 2') 
% %         end
% %         if size(tfdirect,2) ~= 3
% %            error('Test 2: ERROR 3') 
% %         end
% %     end
% %     
% %     
% %     
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Test 3: results should be of size 4 
% % % 
% % % difforder = 1
% % 
% % 
% %     Sinputdata = struct('coordinates',sources);
% %     Sinputdata.doaddsources = 1;
% %     Rinputdata = struct('coordinates',receivers(1,:));
% %     controlparameters.frequencies = frequencies;
% %     controlparameters.difforder = 1;
% %     
% %     EDmain_convexESIE(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);
% % 
% %     eval(['load ',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tf.mat tfdirect tfdiff'])
% % 
% %     if ndims(tfdirect) > 2
% %         error('Test 3: ERROR 1')
% %     else
% %         if size(tfdirect,1) ~= 4
% %            error('Test 3: ERROR 2') 
% %         end
% %         if size(tfdirect,2) ~= 1
% %            error('Test 3: ERROR 3') 
% %         end
% %     end
% %     
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Test 4: results should be of size 4 * 3 
% % % 
% % % difforder = 1
% % 
% % 
% %     Sinputdata = struct('coordinates',sources(1,:));
% %     Sinputdata.doaddsources = 1;
% %     Rinputdata = struct('coordinates',receivers);
% %     controlparameters.frequencies = frequencies;
% %     controlparameters.difforder = 1;
% %     
% %     EDmain_convexESIE(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);
% % 
% %     eval(['load ',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tf.mat tfdirect tfdiff'])
% % 
% %     if ndims(tfdirect) > 2
% %         error('Test 4: ERROR 1')
% %     else
% %         if size(tfdirect,1) ~= 4
% %            error('Test 4: ERROR 2') 
% %         end
% %         if size(tfdirect,2) ~= 3
% %            error('Test 4: ERROR 3') 
% %         end
% %     end
% %     
% %     
% %     




