function [alpha,F] = absolute_ratio(af,H,h)
%ABSOLUTE_RATIO Calculates the absolute ratio and its routing function.
%
% [alpha,F]=absolute_ratio(af[, H[, h]])
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
% af                      - uarcflow object
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

[alpha,F]=oblivous_ratio1(af,H,h);

end
