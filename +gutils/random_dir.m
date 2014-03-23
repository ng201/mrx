function adj = random_dir(n,p)
%RANDOM_DIR Generates a random directed graph.
%
% function adj = random_dir(n,p)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Generates a random directed graph.
%
% Syntax:
%
% adj = gutils.random_dir(10,.2)
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% N                       - number of nodes
% p                       - probability of edge insertion; default is 0.5
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% adj                     - adjacency matrix
%
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

adj=zeros(n); % initialize adjacency matrix

error(nargchk(1,2,nargin));
if nargin==1; 
    p=0.5; 
end;

for i=1:n
  for j=1:i-1
    if rand<=p; adj(i,j)=1; end;
  end
  for j=i+1:n
    if rand<=p; adj(i,j)=1; end;
  end
end

end

