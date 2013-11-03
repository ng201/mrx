classdef mrx
%MRX MRX option object (structure).
%
% function mos=mrx(varargin)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Stores and loads an MRX option structure. The value is stored in a
% persistent variable.
%
% Syntax:
%
% mos = mrx('param1','param2')
% mrx
% mrx('param1','param2')
% mos = mrx
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% varargin                - the arguments in param-value form
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
% mos                     - the mrx options structure
%
%            
% Copyright is with the following author(s):
%
% (C) 2013 Gábor Németh, BME TMIT
%          nemethgab@tmit.bme.hu
% ---------------------------------------------------------------------------
%%
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
%%

        
    properties (SetAccess = 'private')
       
        version = '1.0.0';
        
        OptTol = 1e-8;   % Default absolute tolerance
        EpsTol = 1e-5;   % Default epsilon value
        Infty = 1e9;     % Default infinity value
        Model = 'FULL'; 
        MaxIter = -1;
        TimeLimit = -1;
        
        LPSolver = 0;         % the LP solver ins use
        
        %% Default LP Solver        
        % ------------------
        % which LP solver to use:
        %   0 - 'auto'
        %   1 - 'linprog'   - Matlab's linprog
        %   2 - 'glpk'      - GLPK solver
        %   3 - 'cplexint'  - CPLEX 9 LP solver
        %   4 - 'cplex'     - CPLEX interfaced by IBM 
        %   5 - 'cdd cc'    - CDD Criss-Cross
        %   6 - 'cdd ds'    - CDD Dual Simplex
        %   7 - 'mosek'     - MOSEK
        %
        % 901 - 'mpt2'      - the solver of the Multi Parametric Toolbox v2
        % 902 - 'mpt'       - the solver of the Multi Parametric Toolbox v3
        %%
        
        
        ExtremeSolver = 0;    % the extreme solver in use 
        
        %%        
        % 0 - 'auto'
        % 1 - 'extrpts'
        % 2 - 'lrs'
        % 5 - 'cdd'
        %%
        
        LPFun = [];       % handle to the LP solver
        ExtFun = [];      % handle to the extreme solver         
        
        testSolvers = 1;     % whether test solvers on startup
        
        MsgLevel = 0;
        
        options = struct();        % options passed to the solver
        
    end

    
    methods (Static,Access=private)
        
        out=sI(sclass,stype)
            
    end
    
    methods
        
        function this=mrx(varargin)       
            
            this.options = struct();
            
            % some parameters are removed from input if available 
            showversion = 'off';           
            verbose = 0;
            
            if(nargin > 0)
    
                for i=nargin:-1:1
                    % show version
                    if strcmpi(varargin{i},'version')
                        showversion = 'on';
                        varargin(i) = [];                   
                    elseif strcmpi(varargin{i},'verbose')
                        verbose = 1;
                        varargin(i) = [];
                    end
                end
            end
            
            if verbose
                this.print
            end

            % parse all the other argument with inputParser
            input_parser=inputParser;
            input_parser.CaseSensitive = true;

            addOptional(input_parser,'mos',[]);
            addParamValue(input_parser,'OptTol',1e-8,@isnumeric);
            addParamValue(input_parser,'EpsTol',1e-5,@isnumeric);
            addParamValue(input_parser,'Infty',1e+9,@(x) validateattributes(x, {'numeric'},{'integer','positive'}));
            addParamValue(input_parser,'MaxIter',-1,@(x) validateattributes(x, {'numeric'},{'integer'}));
            addParamValue(input_parser,'TimeLimit',-1,@(x) validateattributes(x, {'numeric'},{'integer'}));

            % select 'auto' to chose the best available
            addParamValue(input_parser,'LPSolver','auto',@(x)any(validatestring(x,{'auto','linprog','glpk','glpkmex','glpkcc','cplex','Cplex','cplexint','cdd','cdd criss-cross','cdd dual-simplex','mpt2','mpt3','mpt'})));
            addParamValue(input_parser,'ExtremeSolver','auto',@(x)any(validatestring(x,{'auto','extrpts','cdd','lrs','.cdd'})));            
            addParamValue(input_parser,'ShowVersion',showversion,@(x)any(validatestring(x,{'on','off'})));
            addParamValue(input_parser,'testSolvers','on',@(x)any(validatestring(x,{'on','off'})));
            addParamValue(input_parser,'MsgLevel','some',@(x)any(validatestring(x,{'on','some','off'})));

            parse(input_parser,varargin{:});           
            
            if isempty(find(strcmp('mos',input_parser.UsingDefaults),1))
                if ~isa(input_parser.Results.mos,'mrx')
                    error('Expected argument %d to be an MRX option object.', i);
                else
                    this=input_parser.Results.mos;
                end
            end   
            
            
            if input_parser.Results.Infty <= input_parser.Results.EpsTol || input_parser.Results.Infty <= input_parser.Results.OptTol
                error('Infty is less or equal to Epstol or OptTol!');
            end                                               
            %%                          
            if ~isempty(find(strcmp('mos',input_parser.UsingDefaults),1)) || isempty(find(strcmp('LPSolver',input_parser.UsingDefaults),1))               
                switch input_parser.Results.LPSolver
                    case 'linprog'
                                                        this.LPSolver=1;
                    case {'glpk','glpkmex','glpkcc'}
                                                        this.LPSolver=2;
                    case 'cplexint'
                                                        this.LPSolver=3;
                    case {'cplex','Cplex'}  
                                                        this.LPSolver=4;
                    case {'cdd','cdd criss-cross'}
                                                        this.LPSolver=5;
                    case 'cdd dual-simplex' 
                                                        this.LPSolver=6;
                    case 'mosek'
                                                        this.LPSolver=7;
                    
                    case 'mpt2'
                                                        this.LPSolver=901;
                    case {'mpt','mpt3'}
                                                        this.LPSolver=902;  
                    otherwise
                                                        this.LPSolver=0;
                end
            end          
            
            %%             
            if ~isempty(find(strcmp('mos',input_parser.UsingDefaults),1)) || isempty(find(strcmp('ExtremeSolver',input_parser.UsingDefaults),1))
                switch input_parser.Results.ExtremeSolver                  
                    case 'cdd'
                                                        this.ExtremeSolver=5;
                    case 'extrpts' 
                                                        this.ExtremeSolver=1;
                    case 'lrs'
                                                        this.ExtremeSolver=2;                                   
                    otherwise
                                                        this.ExtremeSolver=0;
                end
            end
    
            %%
            if ~isempty(find(strcmp('mos',input_parser.UsingDefaults),1)) || isempty(find(strcmp('testSolvers',input_parser.UsingDefaults),1))
                if strcmp(input_parser.Results.testSolvers,'off')
                    this.testSolvers=0;
                else
                    this.testSolvers=1;
                end
            end
                      
            %%
            %---------------------------
            % Default absolute tolerance
            %---------------------------

            % absolute tolerance
            if ~isempty(find(strcmp('mos',input_parser.UsingDefaults),1)) || isempty(find(strcmp('OptTol',input_parser.UsingDefaults),1))
                this.OptTol = input_parser.Results.OptTol;
            end

            %---------------------------
            % Default epsilon value
            %---------------------------

            % epsilon value
            if ~isempty(find(strcmp('mos',input_parser.UsingDefaults),1)) || isempty(find(strcmp('EpsTol',input_parser.UsingDefaults),1))
                this.EpsTol = input_parser.Results.EpsTol;
            end
            %---------------------------
            % Default infinity value
            %---------------------------

            % infinity
            if ~isempty(find(strcmp('mos',input_parser.UsingDefaults),1)) || isempty(find(strcmp('Infty',input_parser.UsingDefaults),1))
                this.Infty = input_parser.Results.Infty;
            end
                        
            %---------------------------
            % Iteration limit
            %---------------------------

            % absolute tolerance
            if ~isempty(find(strcmp('mos',input_parser.UsingDefaults),1)) || isempty(find(strcmp('MaxIter',input_parser.UsingDefaults),1))
                this.MaxIter = input_parser.Results.MaxIter;
            end
            
            %---------------------------
            % Time limit
            %---------------------------

            % absolute tolerance
            if ~isempty(find(strcmp('mos',input_parser.UsingDefaults),1)) || isempty(find(strcmp('TimeLimit',input_parser.UsingDefaults),1))
                this.TimeLimit = input_parser.Results.TimeLimit;
            end              
            
            %%            
            if this.testSolvers || this.LPSolver==0
                % ha nem kell tesztelni, nem ir ki semmit
                [success,solver]=tsL(this,this.testSolvers);
                
                if success==true
                    this.LPSolver=this.sI('lp',solver);
                end                              
                
                if strcmp(solver,'linprog'),
                    disp('WARNING: you have chosen linprog as a default LP solver.');
                    disp('This solver is very slow and numerically not robust!');
                    disp('We strongly advice you to use some other alternative if possible.');
                    fprintf('\n');
                end
                
                switch this.LPSolver
                    case 1
                                switch this.MsgLevel
                                    case 'off'
                                                  displ = 'off';
                                    case 'on'
                                                  displ = 'on';
                                    otherwise
                                                  displ = 'iter';
                                end
                                if this.MaxIter > 0
                                    this.options = optimset('Display',displ,'TolFun',this.OptTol,'MaxIter',this.MaxIter);
                                else
                                    this.options = optimset('Display',displ,'TolFun',this.OptTol);                                    
                                    end                                
                    case 2
                                this.options.tolobj=this.OptTol;
                                this.options.msglev=this.MsgLevel;
                                this.options.itlim=this.MaxIter;
                                this.options.tmlim=this.TimeLimit;
                                
                    case 3
                                this.options.verbose=this.MsgLevel;
                        
                end
                
            end
            
            %%
            if this.testSolvers || this.ExtremeSolver==0
                [success,solver]=tsE(this,this.testSolvers);
                
                if success==true
                    this.ExtremeSolver=this.sI('extreme',solver);
                end 
                
                if strcmp(solver,'extrpts'),
                    disp('WARNING: you have chosen extrpts as a default extreme solver.');
                    disp('This solver is very slow and numerically not robust!');
                    disp('We strongly advice you to use some other alternative if possible.');
                    fprintf('\n');
                end                
                
            end
            
            
            if nargout==0              
                MRXLoad(this);                
            end
            
        end                
        
        function display(this)             
            fprintf('MRX parameters:\n\n');
            disp([['        OptTol'; '        EpsTol'; '         Infty';'   testSolvers'] ...
                  [': ';': ';': ';': '] ...
                  [num2str([this.OptTol; this.EpsTol; this.Infty; this.testSolvers])]]);
            disp(['      LPSolver: ',repmat(' ', 1, 10-length(this.sI('lp', this.LPSolver))),this.sI('lp', this.LPSolver)]);
            disp([' ExtremeSolver: ',repmat(' ', 1, 10-length(this.sI('extreme', this.ExtremeSolver))),this.sI('extreme', this.ExtremeSolver)]);
            if strcmp(this.Model,'SPARSE') 
                disp('  Matrix model:     SPARSE');
            else
                disp('  Matrix model:       FULL');
            end
            fprintf('\n');
        end
        
       
        function net = createNetwork(this,P,u)
            net = pathflow(P,u);
        end
        
    end
end