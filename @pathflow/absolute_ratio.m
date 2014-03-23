function [alpha,F] = absolute_ratio(pf,H,h)
%OBLIVIOUS_RATIO Calculates the absolute ratio and its routing function.
%
% [alpha,F]=absolute_ratio(pf[, H[, h]])
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Calculates the classic absolute ratio, i.e., the max link overload, for
% the path flow object restricted with the Hx<=h inequalities. 
% If H=[] and h=[] it defaults to oblivious ratio.
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% pf                      - pathflow object
% H                       - lhs of Hx<=h
% h                       - rhs of Hx<=h
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% alpha                   - the absolute ratio
% F                       - the routing function
%
% (C) 2012 Gábor Németh, BME-TMIT
%          nemethgab@tmit.bme.hu

narginchk(1,3);

if nargin == 1
    [alpha,F]=oblivious_ratio(pf);
elseif nargin == 2    
    if size(H,2)~=pf.K
        error('Dimension mismatch.');
    end
    
    [alpha,F]=oblivious_ratio2(pf,H);
    
else
    K=pf.K;          % number of source-destination pairs
    rest=size(h,1);  % number of extra inequalities
    
    if size(H,2)~=K || rest~=size(H,1) || size(h,2)>1
        error('Dimension mismatch.');
    end
    
    if rest==1
    %    [alpha,F]=oblivious_ratio3(pf,H,h);
    %elseif isall(h==0)
    %    [alpha,F]=oblivious_ratio2(pf,H);
    else
    %    [alpha,F]=oblivious_ratio4(pf,H,h);
    end
end

end
