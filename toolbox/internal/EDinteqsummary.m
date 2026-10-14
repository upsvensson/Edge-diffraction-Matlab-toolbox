function [textlines,warninglines] = EDinteqsummary(inteqdiagnostics)
% EDinteqsummary makes a text summary of the diagnostics from the solving
% of the integral equation (the struct inteqdiagnostics from
% EDintegralequation_convex_tf), for the log file.
%
% Input parameters:
%   inteqdiagnostics    The struct EDres.inteqdiagnostics
%
% Output parameters:
%   textlines           A cell array with text strings: the solution
%                       method, the number of unknowns and, for the GMRES
%                       methods, the range of the number of iterations and
%                       the largest relative residual; for the iterative
%                       method, how small the last term is, relative to
%                       the first term.
%   warninglines        A cell array with text strings. It is empty, unless
%                       GMRES did not converge for some frequency.
%
% If inteqdiagnostics comes from a results file that was saved by an older
% version of the EDtoolbox, some of the information might be missing. Then
% those lines are left out.
%
% Peter Svensson (peter.svensson@ntnu.no) 12 Oct. 2026
%
% [textlines,warninglines] = EDinteqsummary(inteqdiagnostics);

% 12 Oct. 2026 First version, by Claude, prompted by PS
% 14 Oct. 2026 A line about the _qfinal file (computed and saved, or reused).

textlines = {};
warninglines = {};

dg = inteqdiagnostics;
if ~isstruct(dg) || ~isfield(dg,'solutionmethod') || ~ischar(dg.solutionmethod)
    return
end

textline = ['Solution method: ',dg.solutionmethod];
if isfield(dg,'nQ') && isnumeric(dg.nQ) && isscalar(dg.nQ) && ~isnan(dg.nQ)
    textline = [textline,'. ',int2str(dg.nQ),' unknowns'];
end
if isfield(dg,'nnzH') && isnumeric(dg.nnzH) && isscalar(dg.nnzH) && ~isnan(dg.nnzH)
    textline = [textline,', ',int2str(dg.nnzH),' non-zero matrix elements'];
end
textlines{end+1} = textline;

if strncmp(dg.solutionmethod,'GMRES',5)
    if isfield(dg,'gmresiter') && isfield(dg,'gmresrelres') && isfield(dg,'gmresflag') && ...
            ~isempty(dg.gmresiter) && isequal(size(dg.gmresiter),size(dg.gmresflag))
        niter  = dg.gmresiter(:);
        relres = dg.gmresrelres(:);
        flags  = dg.gmresflag(:);
        ivdone = find(~isnan(flags));
        if ~isempty(ivdone)
            textlines{end+1} = ['GMRES: ',int2str(min(niter(ivdone))),' to ',int2str(max(niter(ivdone))),...
                ' iterations. Largest relative residual: ',num2str(max(relres(ivdone)),3)];
            ivfailed = find(flags ~= 0 & ~isnan(flags));
            if ~isempty(ivfailed)
                [ifreqfailed,~] = ind2sub(size(dg.gmresflag),ivfailed(1));
                warninglines{end+1} = ['WARNING: GMRES did not converge for ',int2str(length(ivfailed)),' of ',...
                    int2str(length(ivdone)),' solutions (the first one: frequency no. ',int2str(ifreqfailed),...
                    '). See EDres.inteqdiagnostics.gmresflag'];
            end
        end
    end
    if isfield(dg,'solvetime') && ~isempty(dg.solvetime) && any(~isnan(dg.solvetime(:)))
        textline = ['GMRES solving time: ',num2str(sum(dg.solvetime(:),'omitnan'),4),' s in total'];
        if isfield(dg,'assemblytime') && ~isempty(dg.assemblytime) && any(~isnan(dg.assemblytime(:)))
            textline = [textline,'. Building the sparse matrix: ',num2str(sum(dg.assemblytime(:),'omitnan'),4),' s in total'];
        end
        textlines{end+1} = textline;
    end
else
    if isfield(dg,'normQ') && isnumeric(dg.normQ) && size(dg.normQ,1) >= 2
        normQ = reshape(dg.normQ,size(dg.normQ,1),[]);
        ivok = find(normQ(1,:) > 0);
        if ~isempty(ivok)
            lasttermratio = max(normQ(end,ivok)./normQ(1,ivok));
            textlines{end+1} = ['Iterative method, ',int2str(size(normQ,1)),' terms: the norm of the last term is at most ',...
                num2str(lasttermratio,3),' of the norm of the first term'];
        end
    end
end

if isfield(dg,'qfinalfile') && ischar(dg.qfinalfile) && ~isempty(dg.qfinalfile)
    textlines{end+1} = ['Integral equation solution (Qfinal): ',dg.qfinalfile];
end
