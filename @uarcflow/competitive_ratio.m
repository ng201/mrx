function [alpha,F] = competitive_ratio(af)
%ABSOLUTE_RATIO Calculates the competitive ratio and its routing function.
%
% [alpha, F]= competitive_ratio(af)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Calculates the clssic oblivous ratio.
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

[alpha,F]=oblivous_ratio1(af);

end