classdef pathflow
% Copyright is with the following author(s):
%
% (C) 2013 Gábor Németh, BME TMIT
%          nemethgab@tmit.bme.hu

% ---------------------------------------------------------------------------
% Legal note:
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
    
    properties (SetAccess = 'private')
        
        P = {};    % paths
        u = [];    % capacity vector
        
        K = 0;     % source-dest pairs
        M = 0;     % edge number
        
        or = 0;    % oblivous ratio if calculated already
        th = {};   % thpol if calculated already
    end
    
    methods
        
        function obj=pathflow(P,u)
            if(nargin > 0)
                if ~iscell(P)
                    error 'pathflow: P must be a cell array';
                end
                if isempty(P)
                    error 'pathflow: P is empty';
                end

                obj.K = length(P);
                obj.M = size(P{1}, 1);

                if size(u, 1) ~= obj.M
                    error 'pathflow: confusing number of links';
                end
            
                obj.P=P;
                obj.u=u;
            end
        end               
        
        function display(pf)
            
            paths=0;
            for i=1:pf.K
                paths=paths+size(pf.P{i},2);
            end
            
            if paths+pf.K<40 & pf.M<=20                
                lb=repmat('[',pf.M,1);
                rb=repmat(']',pf.M,1);
                sp=repmat(' ',pf.M,1);
                ps=repmat('|',pf.M,1);
                xs=sp;
                le=sp;
                xs(floor(pf.M/2)+1)='P';
                le(floor(pf.M/2)+1)='=';
                
                PP=num2str(pf.P{1});
                for i=2:pf.K
                    PP=[PP sp ps sp num2str(pf.P{i})];
                end
                
                disp([xs sp le sp lb sp PP sp rb]);
            else
                disp('P = ');
                disp(pf.P);
            end
            disp(' ');
            if pf.M<20
                disp(['u = [' num2str(pf.u') ']'] );
            else
                disp(['u = [ 1x' num2str(length(pf.u)) ' double ]']);
            end
        end
        
    end
end