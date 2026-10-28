function Q_firstterm = EDinteg_soutermsources(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,controlparameters,...
    Sdata,frequency,ifreq,gaussvectors,sourcelist,addsources,showtext)
% EDinteg_soutermsources - The source term of the edge source integral equation, for one or several sources.
% Calls EDinteg_souterm for each source in sourcelist (for a polygon
% piston: for each of its quadrature points, weighted), and multiplies by
% the source amplitude at the frequency number ifreq. With addsources = 1,
% the contributions of all the sources in sourcelist are added; with
% addsources = 0, sourcelist must be a single source.
%
% Used by EDintegralequation_convex_tf (PART 2) and EDrefineQvec.
%
% Input parameters:
%   envdata, edgedata, edgetoedgedata, Hsubmatrixdata, controlparameters
%                   structs, see help EDdatastructs. Hsubmatrixdata gives
%                   the layout of the source term vector.
%   Sdata           struct, see help EDdatastructs
%   frequency       the frequency, in Hz
%   ifreq           the column of Sdata.sourceamplitudes for this frequency
%   gaussvectors    cell with the Gauss-Legendre nodes and weights for each
%                   number of points (see EDintegralequation_convex_tf)
%   sourcelist      the source numbers
%   addsources      1: the sum over sourcelist. 0: one source.
%   showtext        (optional, default 0)
%
% Output parameters:
%   Q_firstterm     [nQ,1] the source term
%
% See also EDinteg_souterm, EDintegralequation_convex_tf, EDrefineQvec
%
% Peter Svensson (peter.svensson@ntnu.no) 28 Oct. 2026
%
% Q_firstterm = EDinteg_soutermsources(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,controlparameters,...
%    Sdata,frequency,ifreq,gaussvectors,sourcelist,addsources,showtext);

% 28 Oct. 2026 First version (by Claude, prompted by PS): the two copies
% of this code in EDintegralequation_convex_tf (CASE 1 and CASE 2 of
% PART 2), moved here unchanged, so that EDrefineQvec can use it too.

if nargin < 12
    showtext = 0;
end

if addsources == 1
    % CASE 1 of EDintegralequation_convex_tf: the sources are added
    Q_firstterm = 0;
    for isou = sourcelist(:).'
        if strcmp(Sdata.sourcetype,'monopole') == 1
             Q_firstterm_addition = EDinteg_souterm(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,...
                controlparameters,Sdata.vispartedgesfroms(:,isou),Sdata.vispartedgesfroms_start(:,isou),...
                Sdata.vispartedgesfroms_end(:,isou),frequency,gaussvectors,Sdata.rSsho(Sdata.reftoshortlistS(:,isou)),...
                Sdata.thetaSsho(Sdata.reftoshortlistS(:,isou)),Sdata.zSsho(Sdata.reftoshortlistS(:,isou)),showtext);
        elseif strcmp(Sdata.sourcetype,'polygonpiston') == 1
            Q_firstterm_sumpistongaussorder = 0;
            for np = 1:size(Sdata.pistongausscoordinates{1},1)
                Q_firstterm_addition = EDinteg_souterm(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,...
                    controlparameters,Sdata.vispartedgesfroms(:,isou),Sdata.vispartedgesfroms_start(:,isou),...
                    Sdata.vispartedgesfroms_end(:,isou),frequency,gaussvectors,Sdata.rSsho(Sdata.reftoshortlistS(:,isou,np)),...
                    Sdata.thetaSsho(Sdata.reftoshortlistS(:,isou,np)),Sdata.zSsho(Sdata.reftoshortlistS(:,isou,np)),showtext);
                Q_firstterm_sumpistongaussorder = Q_firstterm_sumpistongaussorder + Q_firstterm_addition*squeeze(Sdata.pistongaussweights{isou}(np,1));
            end
            Q_firstterm_addition = Q_firstterm_sumpistongaussorder;
        end
        Q_firstterm = Q_firstterm + Q_firstterm_addition*Sdata.sourceamplitudes(isou,ifreq);
    end
else
    % CASE 2 of EDintegralequation_convex_tf: one source
    if numel(sourcelist) ~= 1
        error('EDinteg_soutermsources:onesource','With addsources = 0, sourcelist must be a single source')
    end
    isou = sourcelist;
    if strcmp(Sdata.sourcetype,'monopole') == 1
        Q_firstterm = EDinteg_souterm(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,...
            controlparameters,Sdata.vispartedgesfroms(:,isou),Sdata.vispartedgesfroms_start(:,isou),...
            Sdata.vispartedgesfroms_end(:,isou),frequency,gaussvectors,Sdata.rSsho(Sdata.reftoshortlistS(:,isou)),...
            Sdata.thetaSsho(Sdata.reftoshortlistS(:,isou)),Sdata.zSsho(Sdata.reftoshortlistS(:,isou)),showtext);
    elseif strcmp(Sdata.sourcetype,'polygonpiston') == 1
        Q_firstterm_sumpistongaussorder = 0;
        for np = 1:size(Sdata.pistongausscoordinates{1},1)
                Q_firstterm_addition = EDinteg_souterm(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,...
                controlparameters,Sdata.vispartedgesfroms(:,isou),Sdata.vispartedgesfroms_start(:,isou),...
                Sdata.vispartedgesfroms_end(:,isou),frequency,gaussvectors,Sdata.rSsho(Sdata.reftoshortlistS(:,isou,np)),...
                Sdata.thetaSsho(Sdata.reftoshortlistS(:,isou,np)),Sdata.zSsho(Sdata.reftoshortlistS(:,isou,np)),showtext);
            Q_firstterm_sumpistongaussorder = Q_firstterm_sumpistongaussorder + Q_firstterm_addition*squeeze(Sdata.pistongaussweights{isou}(np,1));
        end
        Q_firstterm = Q_firstterm_sumpistongaussorder;
    end
    Q_firstterm = Q_firstterm*Sdata.sourceamplitudes(isou,ifreq);
end

end
