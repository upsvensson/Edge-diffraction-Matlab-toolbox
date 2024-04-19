function [q,errbnd] = EDquadgk(integralparameters,a,b,reltol)
% EDquadgk is a stripped-down version of the built-in Matlab function
% quadgk. The integrand, previously in the separate function EDbetaoverml
% is integrated into this function, for efficiency. Optional parameter
% parsing has been removed, and nested functions are mostly included in the
% main code. See "Help quadgk" for more information on the original code.
%
% Input parameters:
%
%   integralparameters  A struct with the fields
%       .rs   .rr   .zs   .zr   .ny     Single values with edge geometry
%       .sinnyfivec  .cosnyfivec        Two arrays of size [1,4] that
%                                       are pi ± theta_S ± theta_R
%       .useterm                        Array of size [1,4] with values 1
%                                       or 0, indicating if the respective
%                                       term should be calculated or not.
%                                       One term might be skipped if there
%                                       is a singularity which is handled
%                                       by the calling function
%                                       (EDwedge1st_fd.m)
%       .Rstart                         The reference distance for the
%                                       phase = 0 definition
%       .k                              The wavenumber
%   a,b                 The start and end points of the integration range.
%   reltol              The relative error tolerance
%
% Output parameters:
%   q                   The integral result
%   errbnd              An estimate of the error bound
%
%   Copyright 2007-2015 The MathWorks, Inc.
%
% Modification by Peter Svensson 19 April 2024 (peter.svensson@ntnu.no)
%
% [q,errbnd] = EDquadgk(integralparameters,a,b,reltol)

% Define constants for the integration
% Comment by PS: In quadgk, the symmetry was employed by mirroring the
% first half but that seemed to take some unnecessary time
%
% Gauss-Kronrod (7,15) pair. Use symmetry in defining nodes and weights.
pnodes = [ ...
    0.2077849550078985; 0.4058451513773972; 0.5860872354676911; ...
    0.7415311855993944; 0.8648644233597691; 0.9491079123427585; ...
    0.9914553711208126];
pwt = [ ...
    0.2044329400752989, 0.1903505780647854, 0.1690047266392679, ...
    0.1406532597155259, 0.1047900103222502, 0.06309209262997855, ...
    0.02293532201052922];
pwt7 = [0,0.3818300505051189,0,0.2797053914892767,0,0.1294849661688697,0];
%NODES = [-pnodes(end:-1:1); 0; pnodes];
NODES = [-0.9914553711208126; -0.9491079123427585; -0.8648644233597691; ...
         -0.7415311855993944; -0.5860872354676911; -0.4058451513773972; ...
         -0.2077849550078985;                   0;  0.2077849550078985; ...
          0.4058451513773972;  0.5860872354676911;  0.7415311855993944; ...
          0.8648644233597691;  0.9491079123427585;  0.9914553711208126];
%WT = [pwt(end:-1:1), 0.2094821410847278, pwt];
WT = [ 0.02293532201052922, 0.06309209262997855, 0.1047900103222502, ...
       0.1406532597155259,  0.1690047266392679,  0.1903505780647854, ...
       0.2044329400752989,  0.2094821410847278,  0.2044329400752989, ...
       0.1903505780647854,  0.1690047266392679,  0.1406532597155259, ...
       0.1047900103222502,  0.06309209262997855, 0.02293532201052922];
%EWT = WT - [pwt7(end:-1:1), 0.4179591836734694, pwt7];
EWT = WT - [0,  0.1294849661688697,  0, 0.2797053914892767, 0, ...
                0.3818300505051189,  0, 0.4179591836734694, 0, ...
                0.3818300505051189,  0, 0.2797053914892767, 0, ...
                0.1294849661688697,  0];
% Fixed parameters.
DEFAULT_DOUBLE_ABSTOL = 1.e-10;
DEFAULT_SINGLE_ABSTOL = 1.e-5;
DEFAULT_DOUBLE_RELTOL = 1.e-6;
DEFAULT_SINGLE_RELTOL = 1.e-4;
MININTERVALCOUNT = 10; % Minimum number subintervals to start.

MAXINTERVALCOUNT = 650;
RTOL = reltol;
ATOL = [];

% Initialize the FIRSTFUNEVAL variable.  Some checks will be done just
% after the first evaluation.
FIRSTFUNEVAL = true;

A = a;
B = b;

tinterval = [-1,1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call a nested function
[tinterval,pathlen] = split(tinterval,MININTERVALCOUNT);
    
% Initialize array of subintervals of [a,b].
subs = [tinterval(1:end-1);tinterval(2:end)];

% Initialize partial sums.
q_ok = 0;
err_ok = 0;

% Initialize main loop
while true
    % SUBS contains subintervals of [a,b] where the integral is not
    % sufficiently accurate. The first row of SUBS holds the left end
    % points and the second row, the corresponding right endpoints.
    midpt = sum(subs)/2;   % midpoints of the subintervals
    halfh = diff(subs)/2;  % half the lengths of the subintervals
    x = NODES*halfh + midpt;
    x = reshape(x,1,[]);   % function f expects a row vector

    % The line below, calling f1, is replaced by the lines below
    % % % % %     [fx,too_close] = f1(x);
    % % % % %     function [y,too_close] = f1(t)
    tt = 0.25*(B-A)*x.*(3 - x.^2) + 0.5*(B+A);

        %The line below, calling evalFun, is replaced by the lines below it
        % % % % %     [fx,too_close] = evalFun(tt);
        % % % % %     function [fx,too_close] = evalFun(x)
        % Evaluate the integrand.
        if FIRSTFUNEVAL
            % Don't check the closeness of the mesh on the first iteration.
            too_close = false;
            % The line below, calling FUN, is replaced by the lines below it
            % fx = FUN(tt);
            % function outval = EDbetaoverml_fd(zvec,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstart,useterm)
                zvec = tt;
                if any(integralparameters.useterm)
                    m = sqrt( (zvec-integralparameters.zs).^2 + integralparameters.rs^2 );
                    l = sqrt( (zvec-integralparameters.zr).^2 + integralparameters.rr^2 );                
                    ml = m.*l;                
                    %------------------------------------------------
                    % We would have liked to use the following code:
                    %       y = (ml + (zvec-zs).*(zvec-zr))./(rs*rr);
                    %       expnyeta = real((real(sqrt(y.^2-1)) + y).^ny);
                    %       coshnyeta = 0.5*( expnyeta + 1./expnyeta);
                    % but to save three vectors, zvec will be re-used for
                    % y and expnyeta and coshnyeta                
                    zvec = (ml + (zvec-integralparameters.zs).*(zvec-integralparameters.zr))./(integralparameters.rs*integralparameters.rr);
                    zvec = real((real(sqrt(zvec.^2-1) + zvec)).^integralparameters.ny);
                    zvec = 0.5*( zvec + 1./zvec);                
                    if any(integralparameters.useterm==0)
                        outval = 0;
                        if integralparameters.useterm(1) == 1
                            outval =  outval + integralparameters.sinnyfivec(1)./(zvec - integralparameters.cosnyfivec(1));
                        end
                        if integralparameters.useterm(2) == 1
                            outval = outval + integralparameters.sinnyfivec(2)./(zvec - integralparameters.cosnyfivec(2));
                        end
                        if integralparameters.useterm(3) == 1
                            outval = outval + integralparameters.sinnyfivec(3)./(zvec - integralparameters.cosnyfivec(3));
                        end
                        if integralparameters.useterm(4) == 1
                            outval = outval + integralparameters.sinnyfivec(4)./(zvec - integralparameters.cosnyfivec(4));         
                        end
                    else
                        outval =          integralparameters.sinnyfivec(1)./(zvec - integralparameters.cosnyfivec(1));
                        outval = outval + integralparameters.sinnyfivec(2)./(zvec - integralparameters.cosnyfivec(2));
                        outval = outval + integralparameters.sinnyfivec(3)./(zvec - integralparameters.cosnyfivec(3));
                        outval = outval + integralparameters.sinnyfivec(4)./(zvec - integralparameters.cosnyfivec(4));                         
                    end                
                    outval = outval./ml.*exp(-1i*integralparameters.k*(m+l-integralparameters.Rstart));
                else
                    outval = 0;
                end            
                fx = outval;
                % END OF INSERTED CALL TO FUN

            finalInputChecks(tt,fx);
            FIRSTFUNEVAL = false;
        else
            too_close = checkSpacing(tt);
            if too_close
                fx = [];
            else
                % The line below, calling FUN, is replaced by the lines below it
                % fx = FUN(tt);
                % function outval = EDbetaoverml_fd(zvec,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstart,useterm)
                    zvec = tt;
                    if any(integralparameters.useterm)
                        m = sqrt( (zvec-integralparameters.zs).^2 + integralparameters.rs^2 );
                        l = sqrt( (zvec-integralparameters.zr).^2 + integralparameters.rr^2 );                    
                        ml = m.*l;                    
                        %------------------------------------------------
                        % We would have liked to use the following code:
                        %       y = (ml + (zvec-zs).*(zvec-zr))./(rs*rr);
                        %       expnyeta = real((real(sqrt(y.^2-1)) + y).^ny);
                        %       coshnyeta = 0.5*( expnyeta + 1./expnyeta);
                        % but to save three vectors, zvec will be re-used for
                        % y and expnyeta and coshnyeta                    
                        zvec = (ml + (zvec-integralparameters.zs).*(zvec-integralparameters.zr))./(integralparameters.rs*integralparameters.rr);
                        zvec = real((real(sqrt(zvec.^2-1) + zvec)).^integralparameters.ny);
                        zvec = 0.5*( zvec + 1./zvec);                    
                        if any(integralparameters.useterm==0)
                            outval = 0;
                            if integralparameters.useterm(1) == 1
                                outval =  outval + integralparameters.sinnyfivec(1)./(zvec - integralparameters.cosnyfivec(1));
                            end
                            if integralparameters.useterm(2) == 1
                                outval = outval + integralparameters.sinnyfivec(2)./(zvec - integralparameters.cosnyfivec(2));
                            end
                            if integralparameters.useterm(3) == 1
                                outval = outval + integralparameters.sinnyfivec(3)./(zvec - integralparameters.cosnyfivec(3));
                            end
                            if integralparameters.useterm(4) == 1
                                outval = outval + integralparameters.sinnyfivec(4)./(zvec - integralparameters.cosnyfivec(4));         
                            end
                        else
                            outval =          integralparameters.sinnyfivec(1)./(zvec - integralparameters.cosnyfivec(1));
                            outval = outval + integralparameters.sinnyfivec(2)./(zvec - integralparameters.cosnyfivec(2));
                            outval = outval + integralparameters.sinnyfivec(3)./(zvec - integralparameters.cosnyfivec(3));
                            outval = outval + integralparameters.sinnyfivec(4)./(zvec - integralparameters.cosnyfivec(4));                         
                        end                    
                        outval = outval./ml.*exp(-1i*integralparameters.k*(m+l-integralparameters.Rstart));
                    else
                        outval = 0;
                    end            
                    fx = outval;
                    % END OF INSERTED CALL TO FUN
            end
        end
        % END OF INSERTED CALL TO evalFun(tt)
    
    if ~too_close
        fx = 0.75*(B-A)*fx.*(1 - x.^2);
    end
    % END OF INSERTED CALL TO f1(x)

    % Quit the while-loop if mesh points are too close.
    if too_close
        break
    end
    fx = reshape(fx,numel(WT),[]);
    % Quantities for subintervals.
    qsubs = (WT*fx) .* halfh;
    errsubs = (EWT*fx) .* halfh;
    % Calculate current values of q and tol.
    q = sum(qsubs) + q_ok;
    tol = max(ATOL,RTOL*abs(q));
    % Locate subintervals where the approximate integrals are
    % sufficiently accurate and use them to update the partial
    % error sum.
    ndx = find(abs(errsubs) <= (2*tol/pathlen)*abs(halfh));
    err_ok = err_ok + sum(errsubs(ndx));
    % Remove errsubs entries for subintervals with accurate
    % approximations.
    errsubs(ndx) = [];
    % The approximate error bound is constructed by adding the
    % approximate error bounds for the subintervals with accurate
    % approximations to the 1-norm of the approximate error bounds
    % for the remaining subintervals.  This guards against
    % excessive cancellation of the errors of the remaining
    % subintervals.
    errbnd = abs(err_ok) + norm(errsubs,1);
    % Check for nonfinites.
    if ~(isfinite(q) && isfinite(errbnd))
        warning(message('MATLAB:quadgk:NonFiniteValue'));
        break
    end
    % Test for convergence; quit the while-loop if fulfilled
    if errbnd <= tol
        break
    end
    % Remove subintervals with accurate approximations.
    % Quit the while-loop if fulfilled
    subs(:,ndx) = [];
    if isempty(subs)
        break
    end
    % Update the partial sum for the integral.
    q_ok = q_ok + sum(qsubs(ndx));
    % Split the remaining subintervals in half. Quit the while-loop if splitting
    % results in too many subintervals. The integral value up to this point
    % will be returned.
    nsubs = 2*size(subs,2);
    if nsubs > MAXINTERVALCOUNT
        warning(message('MATLAB:quadgk:MaxIntervalCountReached',sprintf('%9.1e',errbnd),nsubs));
        break
    end
    midpt(ndx) = []; % Remove unneeded midpoints.
    subs = reshape([subs(1,:); midpt; midpt; subs(2,:)],2,[]);
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Nested functions

    function [x,pathlen] = split(x,minsubs)
        % Split subintervals in the interval vector X so that, to working
        % precision, no subinterval is longer than 1/MINSUBS times the
        % total path length. Removes subintervals of zero length, except
        % that the resulting X will always has at least two elements on
        % return, i.e., if the total path length is zero, X will be
        % collapsed into a single interval of zero length.  Also returns
        % the integration path length.
        absdx = abs(diff(x));
        
        pathlen = x(end) - x(1);
        
        if pathlen > 0
            udelta = minsubs/pathlen;
            nnew = ceil(absdx*udelta) - 1;
            idxnew = find(nnew > 0);
            nnew = nnew(idxnew);
            for j = numel(idxnew):-1:1
                k = idxnew(j);
                nnj = nnew(j);
                % Calculate new points.
                newpts = x(k) + (1:nnj)./(nnj+1)*(x(k+1)-x(k));
                % Insert the new points.
                x = [x(1:k),newpts,x(k+1:end)];
            end
        else
            error('ERROR: pathlen <= 0')
        end
        
        % Remove useless subintervals.
        x(abs(diff(x))==0) = [];
        if isscalar(x)
            % Return at least two elements.
            x = [x(1),x(1)];
        end

    end % split

    function finalInputChecks(x,fx)
        % Do final input validation with sample input and outputs to the
        % integrand function.
        % Check classes.
        if ~(isfloat(x) && isfloat(fx))
            error(message('MATLAB:quadgk:UnsupportedClass'));
        end
        % Check sizes.
        if ~isequal(size(x),size(fx))
            error(message('MATLAB:quadgk:FxNotSameSizeAsX'));
        end
        outcls = superiorfloat(x,fx);
        outdbl = strcmp(outcls,'double');
        % Validate tolerances and apply defaults.
        if isempty(RTOL)
            if outdbl
                RTOL = DEFAULT_DOUBLE_RELTOL;
            else
                RTOL = DEFAULT_SINGLE_RELTOL;
            end
        end
        if isempty(ATOL)
            if outdbl
                ATOL = DEFAULT_DOUBLE_ABSTOL;
            else
                ATOL = DEFAULT_SINGLE_ABSTOL;
            end
        end
        % Make sure that RTOL >= 100*eps(outcls) except when
        % using pure absolute error control (ATOL>0 && RTOL==0).
        if ~(ATOL > 0 && RTOL == 0) && RTOL < 100*eps(outcls)
            RTOL = 100*eps(outcls);
            warning(message('MATLAB:quadgk:increasedRelTol', outcls, sprintf( '%g', RTOL )));
        end
        if outdbl
            % Single RTOL or ATOL should not force any single precision
            % computations.
            RTOL = double(RTOL);
            ATOL = double(ATOL);
        end
    end % finalInputChecks

    function too_close = checkSpacing(x)
        ax = abs(x);
        tcidx = find(abs(diff(x)) <= 100*eps(class(x))*max(ax(1:end-1),ax(2:end)),1);
        too_close = ~isempty(tcidx);
        if too_close
            warning(message('MATLAB:quadgk:MinStepSize',num2str(x(tcidx),6)));
        end
    end % checkSpacing

end % quadgk
