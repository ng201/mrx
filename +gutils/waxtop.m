function [adj_matr,nd_coord]=waxtop(lambda, alpha, beta, domain) 
%WAXTOP creates a random network topology by the method suggested by Waxman in 1988.
%
% function [adj_matr,nd_coord] = waxtop(n,p)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Creates a random network topology by the method suggested by Waxman in 1988.
%
% Syntax:
%
% adj = gutils.waxtop(.2)
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% lambda                  - intensity of the Poisson process
% alpha                   - maximal link probability
% beta                    - parameter to control length of the edges. 
%                           Increased beta yields a larger ratio of long 
%                           edges to short edges
% domain                  - bounds for the region. A 4-dimensional vector in
%                           the form [x_min x_max y_min y_max].
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% adj_matr                - adjacency matrix
% nd_coord                - coordinates of the vertices
%
%            
% Copyright is with the following author(s):
%
% (C) 2001, 2005 R. Gaigalas 
% (C) 2001, 2005 I. Kaj
%          detailed documentation at
%          http://www.math.uu.se/research/telecom/software
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

if (nargin<4) % default parameter values
    lambda = 0.6; % intensity of the Poisson process
    alpha = 0.4; % parameter for the link probability
    beta = 0.1; % parameter for the link probability
    domain = [0 10 0 10]; % bounds for the "geografical" domain
end

xmin = domain(1);
xmax = domain(2);
ymin = domain(3);
ymax = domain(4);
clear domain;

% number of points is Poisson distributed
% with intensity proportional to the area
area = (xmax-xmin)*(ymax-ymin);
npoints = poissrnd(lambda*area);
% npoints = 20;

%given the number of points, nodes are uniformly distributed
nd_coord = rand(npoints, 2);
nd_coord(:, 1) = nd_coord(:, 1)*(xmax-xmin)+xmin;
nd_coord(:, 2) = nd_coord(:, 2)*(ymax-ymin)+ymin;

% create a matrix with all possible distances
x_rep = repmat(nd_coord(:, 1), 1, npoints);
y_rep = repmat(nd_coord(:, 2), 1, npoints);
dist_matr = sparse(triu(((x_rep-x_rep').^2 + ...
    (y_rep-y_rep').^2).^0.5, 1));

% create the matrix of probabilities
prob_matr = alpha*spfun('exp', ...
    -dist_matr./(beta*max(max(dist_matr))));

% generate the adjacency matrix
runi = sprand(dist_matr);
adj_matr = double((runi>0) & (runi < prob_matr));

end
