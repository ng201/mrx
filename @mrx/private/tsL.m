function [success,solver]=tsL(this,test)
%MPT_TSL Tests LP solver or selects best LP solver available.
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%
% Tests LP solver or selects best LP solver available (test & select LP solver).
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% success                     - 1 = 'OK', 0 = 'ERROR'
% solver                      - the code of the solver
%
% Copyright is with the following author(s):
%
% (C) 2013 Gábor Németh, BME-TMIT
%          nemethgab@tmit.bme.hu
%%
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
%%

success=false;

if this.LPSolver~=0    
    switch this.LPSolver
    
        case 1 
                        solvers_list = {...
                                        {'linprog',{'linprog.m'},100},...
                                       };
        case 2 
                        solvers_list = {...
                                        {'glpk',{'glpkcc'},80},...
                                       };
        case 3 
                        solvers_list = {...
                                        {'cplex',{'Cplex'},30},...
                                       };                                        
        case 4 
                        solvers_list = {...
                                        {'cplexint',{'cplexint'},30},...
                                       };                  
        case 5 
                        solvers_list = {...
                                        {'cdd criss-cross',{'cddmex'},70},...
                                       };                  
        case 6 
                        solvers_list = {...
                                        {'cdd dual-simplex',{'cddmex'},71},...
                                       };
        case 7 
                        solvers_list = {...
                                        {'mosek',{'msklpopt'},71},...
                                       };
                                   
        case 901
                        solvers_list = {...
                                        {'mpt2',{'mpt_solveLP'},50},...
                                       };
        case 902
                        solvers_list = {...
                                        {'mpt',{'mpt_solve'},31},...
                                       };
        otherwise
                        solvers_list = {};
    end
   
    found_solvers = SubFindSolvers(solvers_list); 
  
    if ~isempty(found_solvers)
        
        solver = found_solvers{1,1};
        status = SubTestSolver(this,solver);
                
        if status==true
            success=true;           
            return;
        end
        
        fprintf('LP solver %s did not passed test. Will restart search.\n\n',solver_name);      
        
    else
        fprintf('LP solver %s not found. Will restart search.\n\n',solvers_list{1}{1});      
    end
end

solvers_list = {...
    {'linprog',{'linprog.m'},100},...
    {'glpk',{'glpkcc'},80},...
    {'cdd dual-simplex',{'cddmex'},71},...   
    {'cdd criss-cross',{'cddmex'},70},...
    {'mpt2',{'mpt_solveLP'},50}, ...
    {'mpt',{'mpt_solve'},49},...
    {'mosek',{'msklpopt'},30},... 
    {'cplex',{'Cplex'},30},...  
    {'cplexint',{'cplexint'},30},...
    }; 

found_solvers = SubFindSolvers(solvers_list);

%% pick up the solver
LPsolvers = zeros(size(found_solvers,1),1);
disp('MRX tests LP solvers found ...');
disp(' ');
for i=1:size(found_solvers,1)
    
    solver = found_solvers{i,1};   
    LPsolvers(i) = SubTestSolver(this,solver); 
   
end
disp(' ');

LPorder = sortrows(found_solvers(logical(LPsolvers),:),2);

if ~isempty(LPorder)
    success=true;
    solver = LPorder(:,1);
    solver = solver{1};
end

end

%% SUBFUNCTIONS
function found_solvers = SubFindSolvers(solvers_list)
% looks on the path for solvers given in "solvers_list" and shows them on
% the desktop

disp('MRX searches for the available/given LP solver(s) ...');
disp(' ');

found_solvers = [];
for i=1:length(solvers_list)
    solver_name = solvers_list{i}{1};
    solver_file = solvers_list{i}{2};
    solver_pref = solvers_list{i}{3};
    solver_index = [];
    for j=1:length(solver_file)
        % instead of "exist" command, we rely on "which" because class
        % methods might not be visible
        if ~isempty(which(solver_file{j}))
            solver_index = [solver_index j];
        end
    end
    % prepare output to be put on display
    if ~isempty(solver_index)
        str = '';
        for j=1:length(solver_file)
            if j<length(solver_file)
                str = [solver_file{j},', ',str];
            else
                str = [str,solver_file{j}];
            end
        end
        % show found solvers on the screen
        fprintf(' %s %s %s \n',solver_name,...
            repmat('.', 1, 60-length(str)-length(solver_name)), str);
        found_solvers = [found_solvers; {solver_name} {solver_pref}];
    end
end

if ~isempty(found_solvers)
    disp(' ');
end

end 


function status = SubTestSolver(this,solver_name)

status=false;
str='failed';
%mos=mrx(this,'execSolvers','off','LPSolver',solver_name);
org_sol=this.LPSolver;
this.LPSolver=this.sI('lp',solver_name);
try
    %[xopt,fval,exitflag,lambda]=MRXSLscf(this,1,[1;-1],[1;1],[],[],[],[],[],[])
    cmd='[xopt,fval,exitflag,lambda]=MRXSLscf(this,1,[1;-1],[1;1],[],[],[],[],[],[]);';
    T = evalc(cmd);
    if exitflag==1
        status = true;
        str='passed';
    end
end

this.LPSolver=org_sol;

fprintf(' testing %s %s %s \n',solver_name,...
        repmat('.', 1, 52-length(str)-length(solver_name)),str);

end