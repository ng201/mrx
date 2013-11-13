function mflow=max_flow(pf, theta0, H, h)
%MAX_FLOW Calculates max flow
%
% mf=max_flow(pf[, theta0[, H, h]])
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Calculates max flows. If theta0 is set, the amount of traffic specified in 
% theta0 has to be routed in the network. The H x <= h inequalities must be true, also. 
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% pf                      - pathflow object
% theta0                  - min. traffic to be routed; optional
% H                       - the rhs of the Hx<=h inequality; optional
% h                       - the lhs of the Hx<=h inequality; optional
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
% mflow                   - the (row vector of) max flows
%
% (C) 2012 Gábor Németh, BME-TMIT
%          nemethgab@tmit.bme.hu

error(nargchk(1,4,nargin));

if nargin == 1
    theta0=zeros(pf.K,1);
    mflow=max_flow1(pf,theta0);
elseif nargin == 2
    if size(theta0,1)~= pf.K || size(theta0,2)~= 1
        error('theta0 must be a column vector of size %d.',pf.K);
    end
    mflow=max_flow1(pf,theta0);
elseif nargin==4
    if size(theta0,1)~= pf.K || size(theta0,2)~= 1
        error('theta0 must be a column vector of size %d.',pf.K);
    end
    if size(h,2)~= 1
        error('h must be a column vector.');
    end
    if size(H,1)~=length(h) || size(H,2)~=pf.K 
        error('H must be a matrix of rows(h) x %d.',pf.K); 
    end
    
    mflow=max_flow2(pf,theta0,H,h);
else
    error('Invalid number of arguments.');
end

end