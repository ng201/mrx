function [dist,P]=dijkstra(adj,s,target)
%DIJKSTRA Calculates the shortest path between the nodes given.
%
% function [dist,P] = diskstra(adj,s,target)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes the shortest path from s to target.
%
% Syntax:
%
% dist = dijkstra(adj,1,2)
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% adj                     - the adjacency matrix 
% s                       - the source node
% target                  - the target node; if target==[], then dist and P 
%                           include all distances and paths from s
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% dist                    - the length of the shortest path
% P                       - the node list of the path
%
%            
% Copyright (c) 2011, Massachusetts Institute of Technology.
% All rights reserved.
% Gergana Bounova, Oct. 5, 2012
%
% Copyright is with the following author(s):
%
% (C) 2014 Gábor Németh
%          Inter–University Centre for Telecommunications and Informatics
%          Kassai u. 26., Debrecen, Hungary
%          nemethgab@tmit.bme.hu
% ---------------------------------------------------------------------------
%% Legal note:
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
% - Redistributions of source code must retain the above copyright notice,
% this list of
%   conditions and the following disclaimer.
% - Redistributions in binary form must reproduce the above copyright
% notice, this list
%   of conditions and the following disclaimer in the documentation and/or
%   other materials provided with the distribution.
% - Neither the name of the Massachusetts Institute of Technology nor the
% names of its
%   contributors may be used to endorse or promote products derived from
%   this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% ---------------------------------------------------------------------------   

narginchk(3,3);

[n,m]=size(adj);

if n~=m
    error('+gutils/dijkstra: matrix dimension must agree.');
end

% initialize distances ==========================
n=length(adj);            % number of nodes
adjL=adj2adjL(adj);       % list of neighbors

dist=inf(1,n);
dist(s)=0;

previous=[1:n; inf(1,n)]';  % {i: inf}, i=1:n, inf -> not assigned
S=cell(1,n); % shortest path sequence


Q=[1:n]; % all unvisited vertices, entire graph
while length(Q)>0 % while not empty
    % get min dist member among unvisited vertices
    [mindist,min_ind]=min(dist(Q));
    u=Q(min_ind);
    
    % termination condition - save source-u path
    S{u}=[];
    t=u;
    while not(isempty(find(previous(:,1)==t)))  % t in previous.keys():
        % insert u at the beginning of S
        S{u}=[t S{u}];
        t=previous(t,2);
    end
    if length(target)>0 & u==target
        dist=dist(u); P=S{u};
        return
    end            
    
    % =========================================
    Q=setdiff(Q,u);  % remove u from Q
    for v=1:length(adjL{u})   % across all neighbors of u
        v=adjL{u}(v);        
        alt=dist(u)+adj(u,v);
        if alt < dist(v)
            dist(v)=alt;
            previous(v,2)=u;
        end
    end
end

P=S;

end

function L = adj2adjL(adj)

L=cell(length(adj),1);

for i=1:length(adj); 
    L{i}=find(adj(i,:)>0); 
end

end