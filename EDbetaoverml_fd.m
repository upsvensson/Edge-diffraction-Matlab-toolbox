function outval = EDbetaoverml_fd(zvec,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstart,useterm)
% EDbetaoverml_fd - Integrand function which is called for num. int. of ED TF.
%
% Input parameters:
%  zvec,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstart,useterm
% 
% Output parameter:
%  The integral value
%
% Peter Svensson 16 Jan. 2018 (peter.svensson@ntnu.no)
%
% outval = EDbetaoverml_fd(zvec,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstart,useterm);

% 18 July 2010 FUnctioning version
% 1 Nov 2016 Introduced the input parameter useterm = vector of 4 values, 0
% or 1. If the value is zero, then that term (of the four) should be set to
% zero due to a singularity.
% 28 Nov. 2017 Copied to EDtoolbox
% 29 Nov. 2017 Name was ESIE2betaoverml_fd
% 16 Jan. 2018 Changed the use of "useterm". Before, numerical problems
% were encountered for the terms that were not supposed to be computed.

if any(useterm)
    m = sqrt( (zvec-zs).^2 + rs^2 );
    l = sqrt( (zvec-zr).^2 + rr^2 );

    ml = m.*l;

    %------------------------------------------------
    % We would have liked to use the following code:
    %       y = (ml + (zvec-zs).*(zvec-zr))./(rs*rr);
    %       expnyeta = real((real(sqrt(y.^2-1)) + y).^ny);
    %       coshnyeta = 0.5*( expnyeta + 1./expnyeta);
    % but to save three vectors, zvec will be re-used for
    % y and expnyeta and coshnyeta

    zvec = (ml + (zvec-zs).*(zvec-zr))./(rs*rr);
    zvec = real((real(sqrt(zvec.^2-1) + zvec)).^ny);
    zvec = 0.5*( zvec + 1./zvec);

    if any(useterm==0)
        outval = 0;
        if useterm(1) == 1
            outval =  outval + sinnyfivec(1)./(zvec - cosnyfivec(1));
        end
        if useterm(2) == 1
            outval = outval + sinnyfivec(2)./(zvec - cosnyfivec(2));
        end
        if useterm(3) == 1
            outval = outval + sinnyfivec(3)./(zvec - cosnyfivec(3));
        end
        if useterm(4) == 1
            outval = outval + sinnyfivec(4)./(zvec - cosnyfivec(4));         
        end
    else
        outval =          sinnyfivec(1)./(zvec - cosnyfivec(1));
        outval = outval + sinnyfivec(2)./(zvec - cosnyfivec(2));
        outval = outval + sinnyfivec(3)./(zvec - cosnyfivec(3));
        outval = outval + sinnyfivec(4)./(zvec - cosnyfivec(4));         

    end

    outval = outval./ml.*exp(-1i*k*(m+l-Rstart));
else
    outval = 0;
end

    % The gradient instead, for theta_R = pi
    % sinthetaShalf = (cosnyfivec(3)+cosnyfivec(4))/2;

    % % % % costhetaShalf = sinnyfivec(2);
    % % % % % outval = 2*zvec*costhetaShalf.*(2-zvec.^2-costhetaShalf^2)./(zvec.^2 - costhetaShalf^2).^2;
    % % % % outval = 2*(rr./l).^2.*sqrt(1+l/rr).*(2-l/rr);
    % % % % % % % % outval = outval./ml.*exp(-j*k*(m+l-Rstart))/rr;
    % % % % outval = outval./l.*exp(-j*k*(l))/rr;
