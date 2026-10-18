function EDdevguide
% EDdevguide - Checklist for making changes to the EDtoolbox.
% This function contains only documentation. Type
%   help EDdevguide
% or just EDdevguide, to read it. It lists the steps to take for each kind
% of change, so that the version handling, the result recycling (hashes),
% the documentation and the tests stay consistent.
%
% =========================================================================
% A. FOR EVERY CHANGE
% =========================================================================
%  1. Work on a clean state: commit (or stash) earlier changes first, so
%     that "git diff" shows only this change.
%  2. Keep a reference result: before changing anything, run one or two
%     typical cases and save EDres (or tfinteqdiff etc.), so that you can
%     compare afterwards.
%  3. Make the change. In each changed function:
%     - add a line to the local history (after the empty line at the end
%       of the header), starting with the date, e.g.
%         % 18 Oct. 2026 Fixed the ... for ...
%     - update the date in the author line of the header
%         % Peter Svensson (peter.svensson@ntnu.no) 18 Oct. 2026
%     - if the input or output parameters changed, update those sections
%       and the call syntax line.
%  4. Decide: DOES THE CHANGE AFFECT THE COMPUTED RESULTS? Compare with the
%     reference results from step 2. Differences at the level of rounding
%     errors (relative 1e-14 or so) do not count. Then follow B or C.
%  5. Run EDcheckheaders (internal), and EDverify. All tests should pass.
%  6. Add an entry in CHANGELOG.md (in toolbox/), under "## Unreleased"
%     (make that heading at the top if it is not there), in one of the
%     sections
%       **Changes that affect the results**   **Speed**
%       **Input and output**   **Log file**   **Documentation**   **Tests**
%  7. Commit:
%       git add -A          (or git add <the new files>)
%       git status          (check the list)
%       git commit -am "Short description of the change"
%
% =========================================================================
% B. CHANGES THAT DO NOT AFFECT THE RESULTS
%    (speed, log text, documentation, refactoring, tests)
% =========================================================================
%  - No new version number, and no change in EDstageversions.
%  - CHANGELOG.md: an entry under "## Unreleased" (step A6).
%  - Saved stage files remain valid, so nothing is recomputed.
%
% =========================================================================
% C. CHANGES THAT AFFECT THE RESULTS (bug fixes, new methods)
% =========================================================================
%  1. EDstageversions (internal): increase the version number of each
%     changed function in the list "versions", and add a comment, e.g.
%       versions.EDwedge1st_fd = 4;  % 4: <what changed>
%     Every stage that uses the function then gets a new hash, so saved
%     stage files are recomputed instead of reused. (The version number
%     of the toolbox is NOT part of any hash.)
%  2. CHANGELOG.md: an entry under "**Changes that affect the results**",
%     saying which cases are affected and how much.
%  3. If the change fixes a bug or adds a method that can be tested:
%     add a test in EDverify (see F).
%  4. A new version (see G), often right away, so that results can be
%     related to a version.
%
% =========================================================================
% D. A NEW INPUT PARAMETER (or a removed one)
% =========================================================================
%  1. EDcheckinputstructs (internal), local function edparamspec: add a
%     row in the table of the right struct:
%       'name', validator, default, applicability, required, description
%     The description is used by EDinputdoc and the cheatsheet; write it
%     as one or two short sentences.
%     WHICH STRUCT?
%     - Parameters that affect the results: controlparameters (or the
%       struct where it belongs physically). All of controlparameters is
%       part of the run hash (EDinputhash), so EDres is recomputed when it
%       changes.
%     - Settings that do not affect the results (text, saving files,
%       parallel workers): filehandlingparameters. These are not hashed
%       (exception: savediff2result).
%  2. If it affects a stage: add it to the hash of that stage. In the
%     stage function, find the struct EDinputdatastruct (or the hash
%     struct) and add the parameter there. Otherwise a saved stage file
%     computed with another value could be reused!
%  3. EDverify test 17: add a row in the table "changes", and a case in
%     the switch below it, with the stages (stage files) that should get
%     a new hash when the parameter changes.
%  4. If the parameter is stored in a data struct: update EDdatastructs.
%  5. The cheatsheet (see H).
%  REMOVED PARAMETER: delete its row in edparamspec, and add its name to
%  the list removedfields of the struct (in EDcheckinputstructs). Users
%  who still give it then get a warning instead of an error.
%
% =========================================================================
% E. A NEW FUNCTION (or a removed one)
% =========================================================================
%  NEW FUNCTION
%  1. Folder: toolbox/ for functions that users call, toolbox/internal/
%     for all others.
%  2. Header in the standard structure (copy one, e.g. EDinputdoc):
%       function ... = EDname(...)
%       % EDname - One line that says what it does.      <- the H1 line
%       % Description, one or more paragraphs.
%       %
%       % Input parameters:
%       %   ...
%       %
%       % Output parameters:
%       %   ...
%       %
%       % See also ...                                   (optional)
%       %
%       % Peter Svensson (peter.svensson@ntnu.no) <date>
%       %
%       % [out1,out2] = EDname(in1,in2);
%                                                        <- empty line
%       % <date> First version.                          <- history
%     Check with EDcheckheaders.
%  3. If a stage function calls it (directly or indirectly): add it to
%     EDstageversions, both in "versions" (version 1) and in the list
%     "stagesubroutines" of that stage. EDverify test 16 checks these
%     lists against the actual calls (matlab.codetools.
%     requiredFilesAndProducts). This also changes the hash of the stage,
%     so its saved files are recomputed once.
%  4. EDshowcalltree('EDmain_convex') shows the whole call tree, if you
%     want to check where it is used.
%  REMOVED FUNCTION
%  1. Search for calls: grep, or EDshowcalltree.
%  2. Remove it from EDstageversions ("versions" and "stagesubroutines").
%  3. Remove it from git, but keep the file on disk if you want:
%       git rm --cached toolbox/internal/EDname.m
%     and move the file out of the toolbox folders (or to a folder listed
%     in .gitignore). A file that stays on the Matlab path can still be
%     called by mistake.
%  4. CHANGELOG.md: mention it.
%
% =========================================================================
% F. A NEW TEST IN EDverify
% =========================================================================
%  1. Increase ntests, add @edtestNN to the list testfunctions, and write
%     the local function edtestNN (copy a similar test). Each test returns
%     a title, a description (what is computed and the tolerance), result
%     lines, and pass/fail (use edresult).
%  2. Add the test to the list in the header of EDverify, and a history
%     line.
%  3. Check that the test FAILS with the old code and PASSES with the new.
%  4. Tests that use result recycling must not depend on files from earlier
%     EDverify runs (see test 18: delete them first).
%
% =========================================================================
% G. A NEW VERSION
% =========================================================================
%  1. EDgetversion (internal): versionnumber, lastchangedate,
%     lastchangetime, and a history line.
%  2. CHANGELOG.md: rename "## Unreleased" to "## 0.xyz (date)", with a
%     one-paragraph summary (does it change results? which?).
%  3. Run EDverify (the log file name then has the new version number).
%  4. Commit and tag:
%       git commit -am "Version 0.xyz: ..."
%       git tag -a v0.xyz -m "EDtoolbox 0.xyz"
%       git push
%       git push origin v0.xyz
%  Rule of thumb: changes of results -> new version right away (third
%  decimal: 0.613 -> 0.614). Larger changes of structure or method ->
%  second decimal (0.61 -> 0.62).
%
% =========================================================================
% H. THE CHEATSHEET (after a change in edparamspec or ednotes)
% =========================================================================
%  1. In Matlab:
%       EDinputdoc('','<EDtoolbox>/doc/EDinputfields.tex')
%  2. Compile doc/EDtoolbox_cheatsheet.tex (pdflatex), and copy the PDF to
%     the top folder of EDtoolbox (which is in git; doc/ is not).
%  3. EDverify test 19 checks that doc/EDinputfields.tex is up to date.
%  Never edit EDinputfields.tex by hand. The program flow box and the tips
%  are written by hand in EDtoolbox_cheatsheet.tex.
%
% =========================================================================
% I. USEFUL TOOLS
% =========================================================================
%   EDverify          all the tests (about 20); run after every change
%   EDcheckheaders    checks the header structure of all functions
%   EDshowcalltree    the call tree of a function
%   EDstageversions   versions of the functions and the stage dependencies
%   EDinputdoc        all input parameters (from edparamspec)
%   EDdatastructs     the internal data structs
%   EDversionhistory  shows CHANGELOG.md
%
% See also EDverify, EDstageversions, EDcheckinputstructs, EDinputdoc,
% EDdatastructs
%
% Peter Svensson (peter.svensson@ntnu.no) 18 Oct. 2026
%
% EDdevguide

% 18 Oct. 2026 First version (by Claude, prompted by PS).

help(mfilename)
