function EDversionhistory
% This function, EDversionhistory, prints out the EDtoolbox version
% history on the screen.
%
% Peter Svensson (peter.svensson@ntnu.no) 24 Jan 2018

% 24 Jan 2018 First version, just after version 0.101

versionnumbers = cell(1);
versiontext = cell(1);
versiondate = cell(1);

versionnumbers{1} = '0.1';
versiondate{1} = '22 Jan 2018';
versiontext{1} = ['First version of the EDtoolbox. Copied from the ESIE2 toolbox ',...
'with new functions for finding first-order paths (EDfindconvexGApaths) and computing first-order tfs (EDmakefirstordertfs). ',...
'In addition, the analytic integration was implemented in EDwedge1st_fd. This version ',...
'passed EDverify tests 1,3-6 but not test 4 (the corner test).'];

versionnumbers{2} = '0.101';
versiondate{2} = '23 Jan 2018';
versiontext{2} = ['A bug was fixed in EDfindconvexGApaths, which gave the wrong ',...
'amplitude of the direct sound for obscured receivers. This happened only in ',...
'special cases where there happened to be an extra edge or corner hit which was ',...
'unrelated to the obscuring plane. EDverify failed for tests 1,2,4 (did not test 6), so ',...
'some error has sneaked into EDwedge1st_fd.'];

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
