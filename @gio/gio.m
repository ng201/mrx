classdef gio < handle
%GIO Graph io handle object for mrx.
%
% function g=gio(varargin)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Stores and loads a graph io structure.
%
% Syntax:
%
% g = gio('param1','param2')
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% varargin                - 
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% g                       - the gio handle object
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
       
        directed = false;   % whether it is directed or undirected
        model = 'pf';       % the model: pf or af (path-flow, arc-flow)
        
        adj = [];           % adjacency matrix
        cap = [];           % capacity matrix
        edges= [];          % edge list with label and capacity
        
        sdpairs = [];       % source-destination pairs
        P = {};             % paths; each row is a path
        
    end
  
    methods (Static, Access = private)

        % Load file written in .i1.j1.k1 format.
        [adj,cap] = load_i1j1k1(fname)
        
        [adj,cap] = load_lgf(fname,varargin)
        
        % Check whether a sequence of number is graphic, i.e., a graph with this degree sequence exists
        function B = is_graphic(seq)
            if not(isempty(find(seq<=0))) || mod(sum(seq),2)==1
                % there are non-positive degrees or their sum is odd
                B = false; return;
            end
            n=length(seq);
            seq=-sort(-seq);  % sort in decreasing order
            for k=1:n-1
                sum_dk = sum(seq(1:k));
                sum_dk1 = sum(min([k*ones(1,n-k);seq(k+1:n)]));
   
                if sum_dk > k*(k-1) + sum_dk1; B = false; return; end
            end
            B = true;
        end

        % Makes a symmetric one from the input matrix.
        function out = symmetrize(in)        
            % Extract the upper triangular part
            %%%B = sparse(triu(in));
            % Undo
            %%%out = full(B);
            %%%out = out + out' - diag(diag(out));
            
            out = max(in,transpose(in));
        end
                
        function el = generate_el(adj,cap)
            %error(nargchk(1,2,nargin));
            narginchk(1,2);
            
            n=length(adj); % number of nodes
            ed=find(adj>0); % indices of all edges
            el=[];
            count=1;
            
            if nargin==1                
                for e=1:length(ed)
                    [i,j]=ind2sub([n,n],ed(e)); % node indices of edge e  
                    el=[el; i j count adj(i,j)];
                    count=count+1;
                end
            else                
                for e=1:length(ed)
                    [i,j]=ind2sub([n,n],ed(e)); % node indices of edge e  
                    el=[el; i j count cap(i,j)];
                    count=count+1;
                end
            end
        end
        
    end
        
    methods
        
        function this = gio(varargin)
            
            % parse all the other argument with inputParser
            input_parser=inputParser;
            input_parser.CaseSensitive = true;

            addParameter(input_parser,'Directed',true,@(x) validateattributes(x, {'logical'},{}));
            addParameter(input_parser,'Input','',@(x) validateattributes(x, {'char'},{}));
            addParameter(input_parser,'Model','pf',@(x)any(validatestring(x,{'pf','af'})));
            addParameter(input_parser,'FileType','.i1.j1.k1',@(x)any(validatestring(x,{'.i1.j1.k1','lgf'})));

            parse(input_parser,varargin{:});
            
            this.directed = input_parser.Results.Directed;
            this.model = input_parser.Results.Model;
                                                
            if isempty(find(strcmp('Input',input_parser.UsingDefaults),1))
                this.load(varargin);                
            end   
                        
        end

        function display(this)
            fprintf('GIO parameters:\n\n');
            
            v=[];
            
            if this.directed
                v=[v;' true'];
            else
                v=[v;'false'];
            end
            if ~isempty(this.adj)
                v=[v;' true'];
            else
                v=[v;'false'];
            end
            if ~isempty(this.sdpairs)
                v=[v;' true'];
            else
                v=[v;'false'];
            end
            v=[sprintf('   %s',this.model);v];
            
            disp([['       Model'; '    Directed'; '   Generated'; '   Populated'] ...
                  [': ';': ';': ';': '] ...
                  v]);
            fprintf('\n');            
        end
    
        function [] = generate(this,n,varargin)
            input_parser=inputParser;
            input_parser.CaseSensitive = true;
            
            addParameter(input_parser,'Directed',this.directed,@(x) validateattributes(x, {'logical'},{}));

            addParameter(input_parser,'Probability',.5,@(x) validateattributes(x, {'numeric'},{'scalar'}));
            addParameter(input_parser,'Distribution','uniform',@(x)any(validatestring(x,{'uniform','normal','binomial','exponential','sequence'})));

            addParameter(input_parser,'DegreeSequence',[],@(x) validateattributes(x, {'double'},{'vector'}));
            addParameter(input_parser,'MaxEdges',inf,@(x) validateattributes(x, {'uint16'},{'scalar'}));
            
            parse(input_parser,varargin{:});
            
            this.directed = input_parser.Results.Directed;
            
            if this.directed == true
                this.adj = gutils.random_dir(n,input_parser.Results.Probability);                
                this.edges=this.generate_el(this.adj);
                this.cap=[];
            else
                this.adj = gutils.random_undir(n,input_parser.Results.Probability,input_parser.Results.MaxEdges,input_parser.Results.Distribution,input_parser.Results.DegreeSequence);
                this.edges=this.generate_el(this.adj);
                this.cap=[];
            end
            
        end
        
        function net = realize(this)
            
            if isempty(this.adj) 
                err = MException('Network error: please generate or load a new network.');
                throw(err); 
            end
                        
            if strcmp(this.model,'pf')                
                % generate_pf
                if isempty(this.P)
                    err = MException('Network error: please populate the network first.');
                    throw(err); 
                end
                net=generate_pf(this);
            else                
                % generate_af
                net=generate_af(this);
            end
            
        end
        
    end
    
end
