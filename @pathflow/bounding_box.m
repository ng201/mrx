function [A,b]=bounding_box(pf,varargin)
%BOUNDING_BOX Calculates the bouding box
%
% function [A,b]=bounding_box(pf)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Calculates the bounding box og the given pathflow object. max flows. If theta0 
% is set, the amount of traffic specified in theta0 has to be routed in the network. 
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% pf                      - pathflow object
% Options                 - 'PositiveOrthant' sets whether >=0 inequalities should be added                           
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
% A                       - the lhs of the inequalities (Ax <=b) describing 
%                           the bounding box
% b                       - the rhs of the inequalities describing the 
%                           bounding box
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
addRequired(input_parser,'pf');
addParamValue(input_parser,'PositiveOrthant','off',@(x)any(validatestring(x,{'on','off'})));
parse(input_parser,pf, varargin{:});

mf=max_flow(pf);

if strcmp(input_parser.Results.PositiveOrthant,'on')
    A=[eye(pf.K); -eye(pf.K)];
    b=[mf; zeros(pf.K,1)];
else
    A=eye(pf.K);
    b=mf;
end

end
