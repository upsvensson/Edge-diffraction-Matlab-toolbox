function EDversionhistory
% This function, EDversionhistory, prints out the EDtoolbox version
% history on the screen.
%
% Peter Svensson (peter.svensson@ntnu.no) 13 Feb 2018

% 8 Feb 2018 First version, just after version 0.101

versionnumbers = cell(1);
versiontext = cell(1);
versiondate = cell(1);

icounter = 1;
versionnumbers{icounter} = '0.1';
versiondate{icounter} = '22 Jan 2018';
versiontext{icounter} = ['First version of the EDtoolbox. Copied from the ESIE2 toolbox ',...
'with new functions for finding first-order paths (EDfindconvexGApaths) and computing first-order tfs (EDmakefirstordertfs). ',...
'In addition, the analytic integration was implemented in EDwedge1st_fd. This version ',...
'passed EDverify tests 1,3-6 but not test 4 (the corner test).'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.101';
versiondate{icounter} = '23 Jan 2018';
versiontext{icounter} = ['A bug was fixed in EDfindconvexGApaths, which gave the wrong ',...
'amplitude of the direct sound for obscured receivers. This happened only in ',...
'special cases where there happened to be an extra edge or corner hit which was ',...
'unrelated to the obscuring plane. EDverify failed for tests 1,2,4 (did not try test 6), so ',...
'some error has sneaked into EDwedge1st_fd.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.102';
versiondate{icounter} = '24 Jan 2018';
versiontext{icounter} = ['In EDwedge1st_fd, the analytical integration was turned off for now, ',...
'because there seems to be a bug. Using the previous brute-force numerical integration ',...
'instead, EDverify passes tests 1,3-6 again.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.103';
versiondate{icounter} = '25 Jan 2018';
versiontext{icounter} = ['A new version of EDwedge1st_fd, with a much simplified analytical integration, ',...
'which seems to give accurate results. EDverify passes tests 1,3-6.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.104';
versiondate{icounter} = '25 Jan 2018';
versiontext{icounter} = ['A new version of EDsubmatrixstructure, to ensure that there are at least two points per edge.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.105';
versiondate{icounter} = '25 Jan 2018';
versiontext{icounter} = ['A new version of EDfindconvexGApaths. Previously, the direct sound was computed even if .directsound was set to 0.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.106';
versiondate{icounter} = '26 Jan 2018';
versiontext{icounter} = ['A new version of EDmain_convextf. When HOD was not calculated, an ',...
    'empty tfinteqdiff should have been saved, but that was not the case. Fixed now. Also, ',...
    'EDwedge1st_fd was fixed. The case useserialexp2 had not implemented.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.107';
versiondate{icounter} = '26 Jan 2018';
versiontext{icounter} = ['Implemented a new input parameter: Sinputdata.doallSRcombinations. ',...
                  'By setting this to zero, only S1 to R1, S2 to R2, etc will be ',...
                  'computed. Also, changed defaults for saving so that saveeddatafile ',...
                  'and saveSRdatafiles = 1, if not specified. Finally, removed the ',...
                  'addition of "results" to the outputdirectory. Now, saving happens ',...
                  'in exactly the directory you specify. If you dont specify an ',...
                  'outputdirectory, a folder called "results" is generated in the ',...
                  'input folder.'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.108';
versiondate{icounter} = '2 Feb 2018';
versiontext{icounter} = ['Corrected the handling of sourceamplitudes; the frequency- ',...
                  'dependence was not handled before. Developed EDdebug to run through ',...
                  'many combinations of nsource, nreceivers, nfrequencies etc. ',...
                  'Added timingstruct as output data in the tfinteq file. Fixed ',...
                  'small bugs in EDfindconvexGApaths'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.109';
versiondate{icounter} = '9 Feb 2018';
versiontext{icounter} = ['Introduced a new parameter: controlparameters.'...
                        'skipfirstorder (default = 0).'];

icounter = icounter + 1;
versionnumbers{icounter} = '0.110';
versiondate{icounter} = '9 Feb 2018';
versiontext{icounter} = ['Introduced inputdatahash in many functions, which makes it',...
                         'it possible to recycle result files. Uses the function DataHash ',...
                         'from Matlab Central.'];
              
icounter = icounter + 1;
versionnumbers{icounter} = '0.111';
versiondate{icounter} = '13 Feb 2018';
versiontext{icounter} = ['Introduced the function lgwt from Matlab Central, for Gauss-Legendre weights, ',...
                         'in EDdistelements. Replaced the quad_gauss.];

icounter = icounter + 1;
versionnumbers{icounter} = '0.112';
versiondate{icounter} = '14 Feb 2018';
versiontext{icounter} = ['Modified the recycling of results files so that the tf and tfinteq files always contain data, ',...
                         'instead of pointing at another file. Also changed the input struct names to geoinputdata, ',...
                         'Sinputdata, and Rinputdata'];

nversions = size(versiontext,2);
txtwidth = 80;

disp(' ')
disp('******************************************************')
disp('Version history of the EDtoolbox')
for ii = 1:nversions
   disp(' ')
   disp(['Version ',versionnumbers{ii},',   ',versiondate{ii}]) 
   txtstr = versiontext{ii};
   nchars = size(txtstr,2);
   txtstrfinished = 0;
   
   while txtstrfinished == 0
       nchars = size(txtstr,2);
       if nchars <= txtwidth
           disp(txtstr);
           txtstrfinished = 1;
       else
            iv = find(txtstr == 32);
            iv2 = find(iv<=txtwidth);
            nstop = iv(iv2(end));
            disp(txtstr(1:nstop-1));
            txtstr = txtstr(nstop+1:end);
       end
   end
end
