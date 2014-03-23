function [alpha,F] = competitive_ratio(pf)
%OBLIVIOUS_RATIO Calculates the oblivious ratio and the oblivious routing function
%
% [alpha,F]=oblivious_ratio(pf)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Calculates the classic oblivious ratio, i.e., solves the 
%       min alpha : [F_1;...;F_K] = F >= 0, 1^T F_k = 1 k=1,...,K
%                   \forall (i,j) \in E: PFx<=alpha*u
% linear program.
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% pf                      - pathflow object
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% alpha                   - the oblivious ratio
% F                       - the routing function
%
% (C) 2012 Gábor Németh, BME-TMIT
%          nemethgab@tmit.bme.hu

[alpha,F]=oblivious_ration(pf);

end

