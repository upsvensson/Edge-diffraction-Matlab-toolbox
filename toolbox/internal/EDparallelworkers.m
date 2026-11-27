function [nworkersuse,needprobe,parallelnote] = EDparallelworkers(nworkers,tserial,tthreshold)
% EDparallelworkers - The number of parallel workers to use for a parfor
% loop: a parallel pool is started only if the work is large enough.
%
% A parfor loop with nworkers > 0 starts a parallel pool, if none is
% running (with Matlab's default setting). That takes 10-30 s, which can
% be much more than the work in the loop, for instance when the solution
% of the integral equation is read from a _qfinal file and only the
% propagation to a few receivers is computed. Therefore:
%   - if a parallel pool is already running, it is used (nworkersuse =
%     nworkers);
%   - if no pool is running, and the serial computing time of the loop is
%     not known (tserial = []), the caller should first run one iteration
%     serially, to estimate it (needprobe = 1, nworkersuse = 0);
%   - with the estimated serial time of the remaining iterations, tserial,
%     a pool is started (nworkersuse = nworkers) only if tserial is larger
%     than tthreshold (default 30 s); otherwise the loop runs serially.
% To use the workers also for short computations: start the pool before
% the call to EDmain_convex (parpool).
%
% Input parameters:
%   nworkers        the largest number of workers, as given by
%                   filehandlingparameters.parallelworkers (0 = serial)
%   tserial         (optional) [] or the estimated serial time, in s, of
%                   the remaining iterations of the loop
%   tthreshold      (optional) the smallest serial time, in s, for which a
%                   pool is started. Default: 30.
%
% Output parameters:
%   nworkersuse     the number of workers to give to parfor (0 = serial)
%   needprobe       1 if one iteration should be run serially first, to
%                   estimate tserial (no pool running, tserial = [])
%   parallelnote    '' or a note on the decision, for the log file
%
% Peter Svensson (peter.svensson@ntnu.no) 27 Nov. 2026
%
% [nworkersuse,needprobe,parallelnote] = EDparallelworkers(nworkers,tserial,tthreshold);

% 27 Nov. 2026 First version (by Claude, prompted by PS): the parallel
% pool was started also when the propagation was very short (e.g., after
% the solution was read from a _qfinal file).

if nargin < 2
    tserial = [];
end
if nargin < 3 || isempty(tthreshold)
    tthreshold = 30;
end

nworkersuse = 0;
needprobe = 0;
parallelnote = '';

if isempty(nworkers) || nworkers <= 0
    return
end

if edpoolrunning
    nworkersuse = nworkers;
    return
end

if isempty(tserial)
    needprobe = 1;
    return
end

if tserial >= tthreshold
    nworkersuse = nworkers;
    parallelnote = sprintf(['The serial computing time was estimated at %.3g s, so a parallel ',...
        'pool was started.'],tserial);
else
    parallelnote = sprintf(['The serial computing time was estimated at %.3g s (less than %g s), ',...
        'so no parallel pool was started.'],tserial,tthreshold);
end

end

%--------------------------------------------------------------------------
function running = edpoolrunning
% 1 if a parallel pool is running. gcp('nocreate') does not start a pool.
running = false;
try
    p = gcp('nocreate');
    running = ~isempty(p);
catch
    running = false;
end
end
