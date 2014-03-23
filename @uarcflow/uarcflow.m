classdef uarcflow
%UARCFLOW creates an undirected arc-flow model.
%
%
% uaf = same(NA, E, u);
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Creates an undirected arf-flow model from (i) the node-arc incidence matrix (NA),
% (ii) a matrix whose columns specify the source nodes and endpoints of session,
% and (iii) the link capacity vector.
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% NA                      - nodearc incidence matrix
% E                       - src-dst matrix
% u                       - capacities
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% uaf                     - a new arc-flow object
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
    
    properties (SetAccess = 'private')
        
        NA = [];   % node-arc incidence matrix
        E = [];    % a matrix whose columns specify the source nodes
                   % and endpoints of sessions (e_k^i = -1 if i = s_k, e_k^i = 1 if i = d_k,
        % e_k^i = 0 othwewise)
        u = [];    % capacity vector
        
        K = 0;     % source-dest pairs
        N = 0;     % node number
        M = 0;     % edge number
        
    end
    
    methods
        
        function obj=uarcflow(NA,E,u)
            if(nargin > 0)
                if isempty(NA) || isempty(E)
                    error 'uarcflow: Na or E is empty';
                end
                
                obj.N = size(NA, 1);
                obj.M = size(NA, 2);
                
                if size(E, 1) ~= obj.N
                    error 'uarcflow: size mismatch';
                end
                
                if size(u, 1) ~= obj.M
                    error 'uarcflow: confusing number of links';
                end
                
                obj.NA=NA;
                obj.E=E;
                obj.u=u;
            end
        end
        
        function display(af)
            
            disp('NA = ');
            disp(af.NA);
            disp(' ');
            disp('E = ');
            disp(af.E);
            disp(' ');
            if af.M < 20
                disp(['u = [' num2str(af.u') ']'] );
            else
                disp(['u = [ 1x' num2str(length(af.u)) ' double ]']);
            end
        end
        
    end
end