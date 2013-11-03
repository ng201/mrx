function [xopt,obj,exitflag,lambda]=MRXSLscf(this,f,A,b,Aeq,beq,lb,ub,x0,options)
%MRXSLscf LP solver for full matrices.
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%
% LP solver for full matrices.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
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

if any(b<-this.Infty)
    error('MRXSLscf: Upper bound exceeded. Most probably the problem is unbounded.');
end

switch this.LPSolver
    
    %% Matlab's linprog   
    case 1
           
                [xopt,obj,exitflag,~,lambda]=linprog(f(:),A,b,Aeq,beq,lb,ub,[],this.options);

    %% GLPK solver interfaced by glpkmex
    case 2            
                ctype=[repmat('S',size(Aeq,1),1); repmat('U',size(A,1),1)];            
                [xopt,obj,exitflag,lambda_struct]=glpk(f,[Aeq; A],[beq; b],lb,ub,ctype,[],1,this.options);
 
                switch exitflag
                    case 5                  % optimal, feasible                 
                            exitflag = 1;
                    case {3,4,110,213}      % infeasible                  
                            exitflag = -2;
                    case {6,214}            % unbounded                  
                            exitflag = -3;
                    case 207                % maxiter reached                   
                            exitflag = 0;
                    case {1,170}            % problem has no feasible solution                   
                            exitflag = -1;
                    otherwise                  
                            exitflag = -1;
                end
                
                 % assign Lagrange multipliers
                if nargout==4
                    lambda = -lambda_struct.lambda;
                    
                    me=size(Aeq,1);
                    [m,n]=size(A);
                    
                    lambda.ineqlin = -lambda_struct.lambda(me+1:me+m);
                    lambda.eqlin = -lambda_struct.lambda(1:me);
                    if ~isempty(lb)       
                        lambda.lower = -lambda_struct.lambda(me+m+1:me+m+n);
                    else
                        lambda.lower = zeros(n,1);
                    end
                    if ~isempty(ub) && isempty(lb)                      
                        lambda.upper(kept_rows.ub) = -lambda_struct.lambda(me+m+1:me+m+n);
                    elseif ~isempty(ub) && ~isempty(lb)                      
                        lambda.upper = -lambda_struct.lambda(me+m+n+1:me+m+n+n);
                    else
                        lambda.upper = zeros(n,1);
                    end 
                end
            
            

    %% CPLEX interfaced by cplexint
    case 3
                if size(beq,1)~=0
                    INDEQ=[1:1:size(Aeq,1)];                
                    [xopt,obj,~,details]=cplexint([],f',[Aeq; A],[beq; b],INDEQ,[],lb,ub,[],[],this.options);
                else
                    [xopt,obj,~,details]=cplexint([],f',A,b,[],[],lb,ub,[],[],this.options);
                end

                how = lower(details.statstring); 
                if strcmp(how, 'optimal') || strcmp(how, 'optimalrelaxed') || ...
                   strcmp(how, 'optimaltol') || strcmp(how, 'integer optimal solution') || ...
                   strcmp(how, 'optimal with unscaled infeasibilities')
                    exitflag = 1;      
                elseif strcmp(R.how,'unbounded or infeasible')
                    exitflag = -5;
                else
                    exitflag = -1;
                end
                
                % assign Lagrange multipliers
                if nargout==4
                    lambda2 = -details.dual;
                    if ~isempty(lambda2)
                        me=size(Aeq,1);
                        [m,n]=size(A);
                        
                        lambda.ineqlin = lambda2(me+1:me+m);
                        lambda.eqlin = lambda2(1:me);
                        if ~isempty(lb)
                            lambda.lower = lambda2(me+m+1:me+m+n);                             
                        else
                            lambda.lower = zeros(n,1);
                        end
                        if ~isempty(ub) && isempty(lb)                         
                            lambda.upper = lambda2(me+m+1:me+m+n);
                        elseif ~isempty(ub) && ~isempty(lb)
                            lambda.upper = lambda2(me+m+n+1:S.me+m+n+n);
                        else
                            lambda.upper = zeros(n,1);
                        end
                    else
                        lambda = lambda2;
                    end     
                end
           
            
 
    %% Fukuda: CDD - Criss-Cross Method     
    case 5           
                me=size(Aeq,1);
                [m,n]=size(A);
                             
                H.A = [Aeq; A;];
                H.B = [beq; b;];
%                 if ~isempty(lb)
%                     % add lower bounds
%                     H.A = [H.A; -I];
%                     H.B = [H.B; lb];
%                 end
%                 if ~isempty(ub)
%                     % add upper bounds
%                     H.A = [H.A; I];
%                     H.B = [H.B; ub];
%                 end

                ilb = (lb==-Inf) | (lb<=-this.Infty);
                iub = (ub==Inf)  | (ub>=this.Infty);
                % store kept rows
                kept_rows.lb = find(~ilb);
                kept_rows.ub = find(~iub);
                if any(~ilb)
                    % put ones at the positions where there is lb/ub
                    Alb = zeros(nnz(~ilb),n);
                    Alb(:,~ilb) = -eye(nnz(~ilb));
                    H.A = [H.A; Alb];
                    H.B = [H.B; -lb(~ilb)];
                end
                if any(~iub)
                    Aub = zeros(nnz(~iub),n);
                    Aub(:,~iub) = eye(nnz(~iub));
                    H.A = [H.A; Aub];
                    H.B = [H.B; ub(~iub)];
                end 

                H.lin = 1:me;
                if isempty(f)
                    H.obj = zeros(1,n);
                else
                    H.obj = f(:)';
                end
                
                Q = cddmex('solve_lp',H);
                
                obj = Q.objlp;
                xopt = Q.xopt;
           
                if Q.how==1               
                    exitflag = 1;
                elseif Q.how<=5               
                    exitflag = -2;
                elseif Q.how<=7               
                    exitflag = -3;
                else                
                    exitflag = -1;
                end
                
                % assign Lagrange multipliers
                if nargout==4   
                    
                    lambda2 = -Q.lambda;
                    
                    lambda.ineqlin = lambda2(me+1:me+m);
                    lambda.eqlin = lambda2(1:me);
                    if ~isempty(lb)  
                        lambda.lower = zeros(n,1);
                        lambda.lower(kept_rows.lb) = lambda(me+m+1:me+m+numel(kept_rows.lb)); 
                    else
                        lambda.lower = zeros(n,1);
                    end
                    if ~isempty(ub) && isempty(lb)
                        lambda.upper = zeros(n,1);
                        lambda.upper(kept_rows.ub) = lambda(me+m+1:me+m+numel(kept_rows.ub));
                    elseif ~isempty(ub) && ~isempty(lb)
                        lambda.upper = zeros(n,1);
                        lambda.upper(kept_rows.ub) = lambda(me+m+numel(kept_rows.lb)+1:me+m+numel(kept_rows.lb)+numel(kept_rows.ub));
                    else
                        lambda.upper = zeros(n,1);
                    end
                end
 
   
    %% Fukuda: Fukuda: CDD - Dual Simplex Method
    case 6
                me=size(Aeq,1);
                [m,n]=size(A);
                             
                H.A = [Aeq; A;];
                H.B = [beq; b;];
%                 if ~isempty(lb)
%                     % add lower bounds
%                     H.A = [H.A; -I];
%                     H.B = [H.B; lb];
%                 end
%                 if ~isempty(ub)
%                     % add upper bounds
%                     H.A = [H.A; I];
%                     H.B = [H.B; ub];
%                 end

                ilb = (lb==-Inf) | (lb<=-this.Infty);
                iub = (ub==Inf)  | (ub>=this.Infty);
                % store kept rows
                kept_rows.lb = find(~ilb);
                kept_rows.ub = find(~iub);
                if any(~ilb)
                    % put ones at the positions where there is lb/ub
                    Alb = zeros(nnz(~ilb),n);
                    Alb(:,~ilb) = -eye(nnz(~ilb));
                    H.A = [H.A; Alb];
                    H.B = [H.B; -lb(~ilb)];
                end
                if any(~iub)
                    Aub = zeros(nnz(~iub),n);
                    Aub(:,~iub) = eye(nnz(~iub));
                    H.A = [H.A; Aub];
                    H.B = [H.B; ub(~iub)];
                end 

                H.lin = 1:me;
                if isempty(f)
                    H.obj = zeros(1,n);
                else
                    H.obj = f(:)';
                end
                
                Q = cddmex('solve_lp_DS',H)
                
                obj = Q.objlp;
                xopt = Q.xopt;
           
                if Q.how==1               
                    exitflag = 1;
                elseif Q.how<=5               
                    exitflag = -2;
                elseif Q.how<=7               
                    exitflag = -3;
                else                
                    exitflag = -1;
                end
                
                % assign Lagrange multipliers
                if nargout==4   
                    
                    lambda2 = -Q.lambda;
                    
                    lambda.ineqlin = lambda2(me+1:me+m);
                    lambda.eqlin = lambda2(1:me);
                    if ~isempty(lb)  
                        lambda.lower = zeros(n,1);
                        lambda.lower(kept_rows.lb) = lambda(me+m+1:me+m+numel(kept_rows.lb)); 
                    else
                        lambda.lower = zeros(n,1);
                    end
                    if ~isempty(ub) && isempty(lb)
                        lambda.upper = zeros(n,1);
                        lambda.upper(kept_rows.ub) = lambda(me+m+1:me+m+numel(kept_rows.ub));
                    elseif ~isempty(ub) && ~isempty(lb)
                        lambda.upper = zeros(n,1);
                        lambda.upper(kept_rows.ub) = lambda(me+m+numel(kept_rows.lb)+1:me+m+numel(kept_rows.lb)+numel(kept_rows.ub));
                    else
                        lambda.upper = zeros(n,1);
                    end
                end
    
    %% mosek        
    case 7
             a = [Aeq; A];
             
             blc   = [Beq; -inf*ones(size(A,1),1)]; 
             buc   = [Beq; B];           
        
             [res] = msklpopt(f,a,blc,buc,lb,ub); 
             
             if res.solstat == MSK_SOL_STA_OPTIMAL
                 exitflag = 1;
                 xopt = res.sol.bas.xx; 
             else
                 xopt = []; obj = []; lambda=[];
                 exitflag = -1;
             end
        
    %% MPTv2
    case 901
            [xopt,obj,lambda,exitflag,~]=mpt_solveLP(f,A,b,Aeq,beq,x0,[],lb,ub);
         
    %% MPTv3
    case 902
            S.quicklp=true; 
            S.f=f; S.A=A; S.b=b; S.Ae=Aeq; S.be=beq; S.lb=lb; S.ub=ub;
            
            R=mpt_solve(S,1);
  
            xopt=R.xopt; obj=R.obj; lambda=R.lambda; 
            
            exitflag=R.exitflag;
            
    otherwise
        xopt=[]; obj=[]; lambda=[]; exitflag=-1;       
        %error('MRXSLscf: unknown LP solver specified!');
end

end