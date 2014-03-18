function [] = populate(this,K,varargin)
%POPULATE add source-destination pairs and optionaly paths to the graph.
%
% function [] = polupale(this,K,varargin)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Add source-destination pairs and optionaly paths to the graph.
%
% Syntax:
%
% g = gio('param1','param2')
%
% g.populate(4) or populate(g,3,'Method','bimodal','Paths',3)
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% this                    - gio object
% K                       - number of source-destination pairs
% varargin                - other parameters, e.g.,
%                           Method - the methodology of srd-dst selection
%                                    (bimodal by default)
%                           Paths - in path-flow model the number of paths
%                                   per src-dst pair
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% this                    - the modified gio object (note that gio is
%                           handle)
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

n=size(this.adj,1);
m=nnz(this.adj);

if K>n*(n-1)
    err = MException('Parameter value mismatch: k cannot be greater than the size of the underlying graph.');
    throw(err);
end

input_parser=inputParser;
input_parser.CaseSensitive = true;

addParameter(input_parser,'Method','bimodal',@(x)any(validatestring(x,{'bimodal'})));
addParameter(input_parser,'Paths',2,@(x) validateattributes(x, {'int16'},{'scalar'}));

parse(input_parser,varargin{:});

if n<=0 || m<=0
    err = MException('Network error: please generate or load a new network.');
    throw(err);
end

this.sdpairs=[];
this.P={};

switch input_parser.Results.Method
    case 'bimodal'
        order=randperm(n*(n-1),K)
        for k=1:K
            row=floor((order(k)-1)/(n-1))+1;
            col=mod(order(k)-1,(n-1))+1;
            if col>=row col=col+1; end
            this.sdpairs=[this.sdpairs; [row col]];
        end
    otherwise
        error('undefined method');
end

if strcmp(this.model,'af') == 1
    return
end

p=input_parser.Results.Paths;
for k=1:K
    this.P{k}=[];
    adj_tmp=this.adj;
    for i=1:p
        
        [~,R]=gutils.dijkstra(adj_tmp,this.sdpairs(k,1),this.sdpairs(k,2));
        
        path=zeros(m,1);
        for l=2:length(R)
            
            e=find(ismember(this.edges(:,1:2),[R(l-1) R(l)],'rows'),1);
            
            path(this.edges(e,3),1)=1; %%%
            adj_tmp(R(l-1),R(l))=0;
        end
        this.P{k}=[this.P{k},path];
    end
end

end