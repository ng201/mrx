function [] = plot(this,varargin)
%PLOT plots the graph underlying the gio object.
%
% function [] = plot(this,varargin)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Plots the graph underlying the gio object.
%
% Syntax:
%
% g = gio('param1','param2')
%
% plot(g)
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% this                    - gio object
% varargin                - other parameters, e.g.,
%                           OrderBy - degree,betweenness, eigen-centrality,
%                                     module,fiedler-function
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
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
addParameter(input_parser,'OrderBy','none',@(x)any(validatestring(x,{'degree','betweenness','eigen-centrality','fiedler-vector','none','all'})));
addParameter(input_parser,'k',3,@(x)validateattributes(x,{'numeric'},{'scalar'}));
parse(input_parser,varargin{:});
            
adj=this.adj;
n=size(adj,1);

% some scaling
if n < 20
    markersize=3;
elseif n < 100
    markersize=2.5;
elseif n < 250
    markersize=ceil(250/n); % scale for plotting purposes
else 
    markersize=0.5;
end

switch input_parser.Results.OrderBy
    
    case 'degree'
        
        Yd=sort_nodes_by_max_neighbor_degree(adj); % degree centrality
                        
        spy(adj(Yd,Yd),'ks',markersize)
        xlabel('ordered by degree','FontWeight','bold')
        axis([1 n 1 n]);
        set(gca,'YDir','normal')
        axis square

    case 'betweenness'
        
        [betw, Yb] = sort(gutils.node_betweenness(adj)); % node betweenness centrality
        
        spy(adj(Yb,Yb),'ks',markersize)
        xlabel('ordered by betweenness','FontWeight','bold')
        axis([1 n 1 n]);
        set(gca,'YDir','normal')
        axis square

        
    case 'eigen-centrality'
        
        [EC,Yec]=sort(gutils.eigencentrality(adj)); % eigen-centrality
        
        spy(adj(Yec,Yec),'ks',markersize)
        xlabel('ordered by eigen-centrality','FontWeight','bold')
        axis([1 n 1 n]);
        set(gca,'YDir','normal')
        axis square
              
    case 'fiedler-vector'
        
        gutils.ssp(adj,floor(input_parser.Results.k),'plot');
        
    case 'none'
        
        spy(adj,'ks',markersize);
        xlabel('connectivity','FontWeight','bold')
        axis([1 n 1 n]);
        set(gca,'YDir','normal')
        axis square
        
    case 'all'
        
        Yd=sort_nodes_by_max_neighbor_degree(adj); % degree centrality
        [betw, Yb] = sort(gutils.node_betweenness(adj)); % node betweenness centrality
        [EC,Yec]=sort(gutils.eigencentrality(adj)); % eigen-centrality
        
        set(gcf,'Color',[1 1 1])
        subplot(2,2,1)
        spy(adj(Yd,Yd),'ks',markersize)
        xlabel('Ordered by degree','FontWeight','bold')
        axis([1 n 1 n]);
        set(gca,'YDir','normal')
        axis square
        
        subplot(2,2,2)
        spy(adj(Yb,Yb),'ks',markersize)
        xlabel('Ordered by betweenness','FontWeight','bold')
        axis([1 n 1 n]);
        set(gca,'YDir','normal')
        axis square
        
        subplot(2,2,3)
        spy(adj(Yec,Yec),'ks',markersize)
        xlabel('Ordered by eigen-centrality','FontWeight','bold')
        axis([1 n 1 n]);
        set(gca,'YDir','normal')
        axis square
        
        subplot(2,2,4)
        spy(adj,'ks',markersize);
        xlabel('Connectivity','FontWeight','bold')
        axis([1 n 1 n]);
        set(gca,'YDir','normal')
        axis square
        
    otherwise
        
        error('What happended?');
        
end

end


function I = sort_nodes_by_max_neighbor_degree(adj)

[deg,~,~]=gutils.degrees(adj); % compute all degrees, use "deg" assuming symmetry
degmat=zeros(size(adj,1),2);  % a nx2 matrix

for x=1:size(adj,1)   % across all nodes
    degmat(x,1)=deg(x);
    if deg(x)==0
        degmat(x,2)=0; 
        continue; 
    end
    nei_inds=kneighbors(adj,x,1);
    if isempty(nei_inds)
        continue;
    else
        degmat(x,2)=max(deg(nei_inds));
    end
    
end

[sortmat,I]=sortrows(degmat);

end

function kneigh = kneighbors(adj,ind,k)

adjk = adj;
for i=1:k-1
    adjk = adjk*adj; 
end;
kneigh = find(adjk(ind,:)>0);

end
