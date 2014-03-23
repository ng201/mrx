function [a,Mf,theta,mflow]=bounding_simplex(pf, varargin)
%BOUNDING_SIMPLEX returns the equation of the bounding simplex constrined 
%                 by the minimum demand theta0 and cutting plane given by Hx<=h.
%
% [M,a,theta]=bounding_simplex(pf[, theta0 [, H, h [,options]]])
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Returns the equation of the bounding simplex in the following form:
%                         ax<=Mf                                       (1)
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% pf                      - pathflow object
% theta0                  - the base point; optional
% Options                 - type 'TryCenter' to place base point in
%                           the middle ('on','off'); default is 'off';
%                           type 'PositiveOrthant' to add >=0 inequalities;
%                           default is off
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% a                       - see equation (1)
% Mf                      - see equation (1)
% theta                   - demand generating exactly Mf
% mflow                   - the max flows relative to theta0 and cutting
%                           planes Hx<=h
%
%
% Copyright is with the following author(s):
%
% (C) 2014 Gábor Németh
%          Inter–University Centre for Telecommunications and Informatics
%          Kassai u. 26., Debrecen, Hungary
%          nemethgab@tmit.bme.hu
% ---------------------------------------------------------------------------
%% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
% 
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% ---------------------------------------------------------------------------   
%%

input_parser=inputParser;
input_parser.CaseSensitive = true;
addRequired(input_parser,'pf');
addOptional(input_parser,'theta0', [] ,@(x) validateattributes(x, {'numeric'}, {'column','nrows',pf.K}));

addOptional(input_parser,'H',[],@(x) validateattributes(x, {'numeric'}, {'2d','ncols',pf.K}));
addOptional(input_parser,'h',[],@(x) validateattributes(x, {'numeric'}, {'column'}));

addParamValue(input_parser,'TryCenter','off',@(x)any(validatestring(x,{'on','off'})));
addParamValue(input_parser,'PositiveOrthant','off',@(x)any(validatestring(x,{'on','off'})));
parse(input_parser,pf, varargin{:});

theta0=input_parser.Results.theta0;

if isempty(theta0)
    theta0=zeros(pf.K,1);
end

H=input_parser.Results.H;
h=input_parser.Results.h;

if isempty(H) && isempty(h)
    [a,Mf,theta,mflow]=bounding_simplex1(pf,theta0,'TryCenter',input_parser.Results.TryCenter,'PositiveOrthant',input_parser.Results.PositiveOrthant);
else
    [a,Mf,theta,mflow]=bounding_simplex2(pf,theta0,H,h,'TryCenter',input_parser.Results.TryCenter,'PositiveOrthant',input_parser.Results.PositiveOrthant);
end

end
