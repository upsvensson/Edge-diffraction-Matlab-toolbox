function outval = EDbetaoverm_fd(zvec,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,m0,useterm)
% EDbetaoverm_fd - Integrand function which is called for num. int. 
% by EDintegratebetaoverm.
%
% Input parameters:
%  zvec,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,m0,useterm
% 
% Output parameter:
%  The integral value
%
% Peter Svensson 28 May 2018 (peter.svensson@ntnu.no)
%
% outval = EDbetaoverm_fd(zvec,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,m0,useterm);

% 26 May 2018 First version. Copied EDbetaoverml_fd and removed the l-part.
% 28 May 2018 Some correction: R0 was still used but should be m0.

if any(useterm)
    m = sqrt( (zvec-zs).^2 + rs^2 );
%     l = sqrt( (zvec-zr).^2 + rr^2 );
% 
%     ml = m.*l;

    ml = m.*sqrt( (zvec-zr).^2 + rr^2 );
    
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

    outval = outval./m.*exp(-1i*k*(m-m0));
else
    outval = 0;
end

