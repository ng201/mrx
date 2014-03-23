function [alpha,F] = oblivious_ratio(af)
%ABSOLUTE_RATIO Calculates the oblivious ratio and its routing function.
%
% [alpha, F]= oblivious_ratio(af)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Calculates the clssic oblivous ratio, i.e., the max link overload, for
% the path flow object restricted with the Hx<=h inequalities. 
% If H=[] and h=[] it defaults to oblivious ratio.
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% af                      - uarcflow object
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% alpha                   - the absolute ratio
% F                       - the routing function
%
% (C) 2012 Gábor Németh, BME-TMIT
%          nemethgab@tmit.bme.hu

narginchk(1,3);

[alpha,F]=oblivous_ratio1(af);

end
