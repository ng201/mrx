function pf = generate_pf(this)
%GENERATE_PF generates a pathflow object from teh data stored in gio object.
%
% function pf = generate_pf(this)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Designes a pathflow object from teh data stored in gio object.
%
% Syntax:
%
% pf = generate_pf(g);
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% this                    - the gio object
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% pf                      - the generated and condifured pathflow object
%            
% Copyright is with the following author(s):
%
% (C) 2014 G�bor N�meth
%          Inter�University Centre for Telecommunications and Informatics
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
    pf=pathflow(this.P,this.edges(:,4));
end