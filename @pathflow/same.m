function tribool=same(pf, H, h, varargin)
%SAME decides whether the polytope given by Hx<h, x>=0 is the same as the
%thpol of pf.
%
% tribool=same(pf, H, h, varargin)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Decides whether the polytope given by Hx<h, x>=0 is the same as the thpol of pf.
% This approach enumerates the extreme points of the polytope Hx<=h, x>=0 and 
% The extreme points of the deal representation of the THPol are
% enumerated.
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% pf                      - pathflow object
% H                       - the lhs of the Hx<=h equation
% h                       - the rhs
% varargin                - Method: primal, dual, successive, gallo
%                         - vertices
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% tribool                 - yes/no/error
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
addOptional(input_parser,'H',[],@(x) validateattributes(x, {'numeric'}, {'2d','ncols',pf.K}));
addOptional(input_parser,'h',[],@(x) validateattributes(x, {'numeric'}));

addParamValue(input_parser,'Method','dual',@(x)any(validatestring(x,{'primal','dual','successive','gallo','gallo-right','gallo-left'})));
parse(input_parser,pf, varargin{:});

if size(H,1)~=size(h,1)
   error('The number of rows of H have to be the same as the number of rows in h.'); 
end

if strcmp(input_parser.Results.Method,'primal')
    tribool=same1(pf,H,h);
elseif strcmp(input_parser.Results.Method,'dual') 
    tribool=same2(pf,H,h);
elseif strcmp(input_parser.Results.Method,'successive') 
    tribool=same3(pf,H,h);
elseif strcmp(input_parser.Results.Method,'gallo')  || strcmp(input_parser.Results.Method,'gallo-right')
    tribool=same4(pf,H,h);
elseif strcmp(input_parser.Results.Method,'gallo-left') 
    tribool=same4(pf,H,h);
else    
    warning('Undefined method. Fallback to dual.');
    tribool=same2(pf,H,h);
end

end