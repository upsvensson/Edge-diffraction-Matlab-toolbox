% EDdebug runs EDtoolbox for several defined tests, to check its handling of
% different numbers of sources, receivers, and frequencies
% as well as values of difforder = 0,1,2, doaddsources = 0,1
% sourceamplitudes = 1 or ones(nfrequencies,nsources)
% 
% Peter Svensson 21 Nov. 2023 (peter.svensson@ntnu.no)

% 31 Jan 2018 First version
% 3 June 2020 Fixed it so that the path folders could have a space in the name.
% 21 Nov. 2023 Runs EDmain_convex instead of EDmain_xxxxxx

[EDversionnumber,changedate,changetime] = EDgetversion;

% Two sources, three receivers, four frequencies
% doaddsources = 0 or 1

% sources = [0.04 0 0.00001;0 0 -0.32001];
sources = [0.04 0 0.00001;0 0 0.0001;-0.04 0 0.0001];%0 0 -0.32001];
receivers = [0 0 -2; 1 1 1];
frequencies = [1 10 20 40];

sources = EDcirclepoints(0.1,2);
sources(:,3) = 0.00001;

% sources = sources(1:2,:);

mfile = mfilename('fullpath');
[infilepath,filestem] = fileparts(mfile);

corners = [     -0.2000   -0.4400   -0.3200
    0.2000   -0.4400   -0.3200
    0.2000    0.2000   -0.3200
   -0.2000    0.2000   -0.3200
   -0.2000   -0.4400         0
    0.2000   -0.4400         0
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
controlparameters.docalctf = 1;
controlparameters.docalctf_ESIEBEM = 0;
controlparameters.docalcir = 0;
controlparameters.Rstart = 0;
envdata.cair = 344;
Sinputdata = struct;

souvar = [1 size(sources,1)];
recvar = [1 size(receivers,1)];
freqvar = [1 length(frequencies)];
addsouvar = [0 1];
diffordervar = [0 1 2];
souamp = [1 2];

casecounter = 0;
for ii = 1:length(souvar)
    Sinputdata.coordinates = sources(1:souvar(ii),:);
    nsources = size(Sinputdata.coordinates,1);
    for jj = 1:length(recvar)
        Rinputdata.coordinates = receivers(1:recvar(jj),:);
        for kk = 1:length(freqvar)
            controlparameters.frequencies = frequencies(1:freqvar(kk));
            nfrequencies = length(controlparameters.frequencies);
            for ll = 1:length(diffordervar)
                controlparameters.difforder = diffordervar(ll);
                for mm = 1:length(addsouvar)
                    Sinputdata.doaddsources = addsouvar(mm);   
                    for nn = 1:length(souamp)
                        if nn == 1
                            Sinputdata.sourceamplitudes = 1;
                        else
                           Sinputdata.sourceamplitudes = ones(nsources,nfrequencies); 
                        end
                        
                        casecounter = casecounter + 1;
                        disp(['Case no. ',int2str(casecounter),': ',int2str(ii),' ',int2str(jj),' ',int2str(kk),' ',int2str(ll),' ',int2str(mm),' ',int2str(nn)])
%                        EDmain_convexESIE(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);
                        EDmain_convex(geoinputdata,Sinputdata,Rinputdata,struct,controlparameters,filehandlingparameters);
                         eval(['load ''',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tf.mat''',' tfdirect tfdiff'])
                        if controlparameters.difforder >= 2
                             eval(['load ''',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tfinteq.mat''',' tfinteqdiff'])
                        else
                            tfinteqdiff = [];
                        end

                        if any(any(any(tfdiff))) && controlparameters.difforder == 0
                            error('ERROR: difforder was set to zero but tfdiff got some result')
                        end
                        if ~any(any(any(tfdiff))) && controlparameters.difforder > 0
                            disp('WARNING: difforder was set > 0 but tfdiff got no result. Check if this is as expected')
                        end
                        if any(any(any(tfinteqdiff))) && controlparameters.difforder <2
                            error('ERROR: difforder was set <2 but tfinteqdiff got some result')
                        end
                        if ~any(any(any(tfinteqdiff))) && controlparameters.difforder > 1
                            error('ERROR: difforder was set > 1 but tfinteqdiff got no result')
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




