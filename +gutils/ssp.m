% Uses the fiedler vector to assign nodes to groups
% INPUTS: adj - adjancency matrix, k - desired number of nodes in groups [n1, n2, ..], [optional]
% OUTPUTs: modules - [k] partitioned groups of nodes
% Other functions used: fiedler_vector.m

function modules = ssp(adj,k,mode)
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

narginchk(2,3);

% find the Fiedler vector: eigenvector corresponding to the second smallest eigenvalue of the Laplacian matrix
fv = gutils.fiedler_vector(adj);
[~,I]=sort(fv);

% depending on k, partition the nodes
if nargin==1
    
    modules{1}=[]; modules{2}=[];
    % choose 2 groups based on signs of fv components
    for v=1:length(fv)
        if fv(v)>0; modules{2} = [modules{2}, v]; end
        if fv(v)<=0; modules{1} = [modules{1}, v]; end
    end
end

if nargin==2

  k = [0 k];
    
  for kk=1:length(k)
        
    modules{kk}=[];
    for x=1:k(kk); modules{kk} = [modules{kk} I(x+k(kk-1))]; end
         
  end

  modules = modules(2:length(modules));
end

if nargin==3 & strcmp(mode,'plot')
    set(gcf,'Color',[1 1 1])
    subplot(1,2,1)
    plot(fv(I),'k.');
    xlabel('index i')
    ylabel('fv(i)')
    title('Sorted fiedler vector')
    axis('tight')
    axis('square')

    subplot(1,2,2)
    spy(adj(I,I),'k')
    title('Sorted adjacency matrix')
else
    warning('@gio/ssp: invalid argument.');
end

end
