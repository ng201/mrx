% The eigenvalues of the Laplacian of the graph
% INPUTs: adjacency matrix
% OUTPUTs: laplacian eigenvalues, sorted

function s=graph_spectrum(adj)
%GRAPH_SPECTRUM Calculates the eigenvalues of the Laplacian of the graph.
%
% function s = grapth_spectrum(adj)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes the eigenvalues of the Laplacian of the graph.
%
% Syntax:
%
% s = graph_spectrum(adj)
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% adj                     - the adjacency matrix 
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% fv                      - laplacian eigenvalues, sorted
%
%            
% Copyright (c) 2011, Massachusetts Institute of Technology.
% All rights reserved.
% GB: Last updated: Oct. 10, 2012
%
% Copyright is with the following author(s):
%
% (C) 2014 Gábor Németh, BME TMIT
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

error(nargchk(1,1,nargin));

[n,m]=size(adj);

if n~=m
    error('+gutils/graph_spectrum: matrix dimension must agree.');
end

[v,D]=eig(gutils.laplacian_matrix(adj));
s=-sort(-diag(D)); % sort in decreasing order

end