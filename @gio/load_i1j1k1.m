function [adj,cap] = load_i1j1k1(fname)
%LOAD_I1J1K1 loads network topology from file written in i1j1k2 format.
%
% function [adj,cap] = load_i1j1k1(fname)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Loads network topology from file written in i1j1k2 format.
%
% Syntax:
%
% [adj,cap] = load_i1j1k1('textgraph.txt');
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% fname                   - file to read
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% adj                     - the adjacency matrix (with edge lengths)
% cap                     - the matrix of capacities (same as adj)
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

fid = fopen(fname);

tline = fgetl(fid);
while ischar(tline)
    
    if strncmp(tline,'#',1) ~= 1
       
        tline = strtrim(tline); % trim
    
        commas=find(tline==',');
        n=length(commas); % number of nodes
        adj=zeros(n); % initialize adjacency matrix

        if commas(1)>1
            % Extract the neighbors of the first node only
            neigh=str(1:commas(1)-1);
            dots=find(neigh=='.');
            for d=1:length(dots)-1; adj(1,str2num(neigh(dots(d)+1:dots(d+1)-1)))=1; end
            adj(1,str2num(neigh(dots(length(dots))+1:length(neigh))))=1;
        end

        % Extract the neighbors of the remaining 2:n nodes
        for i=2:n
            neigh=str(commas(i-1)+1:commas(i)-1);
            if isempty(neigh); continue; end
    
            dots=find(neigh=='.');
            for d=1:length(dots)-1; adj(i,str2num(neigh(dots(d)+1:dots(d+1)-1)))=1; end
            adj(i,str2num(neigh(dots(length(dots))+1:length(neigh))))=1; 
        end
        
        break;    
    end
    
    tline = fgetl(fid);     % next line
end

fclose(fid);

% undirected - make it symmetric
if this.directed == false
    adj=symmetrize(adj);
end

if nargout=2
    cap=adj;
else
    cap=[];
end

end