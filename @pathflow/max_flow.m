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

narginchk(1,4);

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