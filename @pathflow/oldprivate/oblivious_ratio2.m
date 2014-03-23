function [alpha,F] = oblivious_ratio2(pf, H, varargin)
%OBLIVIOUS_RATIO Calculates the oblivious ratio and the oblivious routing function
%                w.r.t. the given ineq Hx<=0
%
% [alpha,F]=oblivious_ratio(pf, varargin)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Calculates the classic oblivious ratio, i.e., solves the 
%       min alpha : [F_1;...;F_K] = F >= 0, 1^T F_k = 1 k=1,...,K
%                   \forall (i,j) \in E: PFx<=alpha*u; Hx<=0
% linear program.
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% pf                      - pathflow object
% Options                 - type 'TolFun' to set tolerance for linprog;
%                           default is 1e-8; set 'Solver' to glpk to use
%                           glpk, cplexint to use cplex through the cplexint 
%                           interface; default in linprog
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% alpha                   - the oblivious ratio
% F                       - the routing function
%
% (C) 2012 Gábor Németh, BME-TMIT
%          nemethgab@tmit.bme.hu

input_parser=inputParser;
input_parser.CaseSensitive = true;
addRequired(input_parser,'pf');
addOptional(input_parser,'H',[],@(x) validateattributes(x, {'numeric'}, {'2d','ncols',pf.K}));

addParamValue(input_parser,'TolFun',1e-8,@isnumeric);
addParamValue(input_parser,'Solver','linprog',@(x)any(validatestring(x,{'linprog','glpk','cplexint'})));
parse(input_parser,pf, varargin{:});

K = pf.K;   % source-dest pairs
M = pf.M;   % number of edges
P = pf.P;   % paths
u = pf.u;   % capacities

Q = size(H,1); % number of rows of H

p=0;                 % number of path
p_k=zeros(1,K);      % pk(k) is the number of path of k-th sorce-dest

for k=1:K 
    p_k(k)=size(P{k},2);
    p=p+p_k(k);
end

% variables:
%    oblivious ratio, f_1_1,...,f_1_p1,...,f_K_1,...f_K_pK, 
%    \forall (i,j) \in E: w^ij_1,...,w^ij_M
%    \forall (i,j) \in E, \forall k \in K: b^ij_k 
%    \forall (i,j) \in E, \forall Q \in K: la^ij_q 

cols = 1 + p + M*M + K*M + M*Q; % all the cols

% col offsets
cols_w=1 + p;               
cols_b=1 + p + M*M;         
cols_la=1 + p + M*M + K*M;  

Aeq=sparse(K,cols);
beq=sparse(K,1);

% \forall k \in K: \sum_1^pk f_k^l = 1
pi=1; 
for k=1:K    
    for l=1:p_k(k)
        Aeq(k,pi+l)=1; 
    end    
    beq(k,1)=1;    
    pi=pi+p_k(k);    
end

rows = M + M * p + M * K;
A=sparse(rows,cols);
b=sparse(rows,1);

% \forall (i,j) \in E: w^ij u <= oblivious ratio
% M sor
for e=1:M
    A(e,1)=-1;
    A(e,cols_w+(e-1)*M+1:cols_w+e*M)=u';
end

% \forall (i,j) \in E, forall k \in K: 0 >=  -b -w^ij Pk 
pi=1;
for e=1:M
    for k=1:K
        for l=1:p_k(k)
            A(M+pi+l-1,cols_w+(e-1)*M+1:cols_w+e*M)=-P{k}(:,l)';
            A(M+pi+l-1,cols_b+(e-1)*K+k)=-1;
        end
        pi=pi+p_k(k);
    end
end

% forall (i,j) \in E, forall k \in K: b +P^ij_k f_k/u_ij <=0
pi=1;
for k=1:K
    for e=1:M
        A(M+M*p+(e-1)*K+k,cols_b+(e-1)*K+k)=u(e); % b        
        A(M+M*p+(e-1)*K+k,pi+1:pi+p_k(k))=P{k}(e,:); % f
        
        A(M+M*p+(e-1)*K+k,cols_la+(e-1)*Q+1:cols_la+e*Q)=-u(e)*H(:,k); % la 
    end
    pi=pi+p_k(k);
end

if strcmp(input_parser.Results.Solver,'glpk')    
    ctype=[repmat('S',size(Aeq,1),1); repmat('U',size(A,1),1)];
    param.tolobj=input_parser.Results.TolFun;
    param.msglev=2; 
    disp(' ');
    [x,alpha,exitflag]=solvers.glpkmex.glpk([1; zeros(cols-1,1)],[Aeq; A],[beq; b],[zeros(cols_b,1); -inf*ones(K*M,1); zeros(M*Q,1)],[],ctype,[],1,param);
    if exitflag~=5
        error('GLPK exited with error code %d.',exitflag); 
    end
elseif strcmp(input_parser.Results.Solver,'cplexint') 
    INDEQ=[1:1:size(Aeq,1)];
    OPTIONS.verbose=1;
    [x,alpha,exitflag,details]=solvers.cplexint.cplexint([],[1; zeros(cols-1,1)],[Aeq; A],[beq; b],INDEQ,[],[zeros(cols_b,1); -inf*ones(K*M,1); zeros(M*Q,1)],[],[],[],OPTIONS);
    if exitflag~=1
        error('CPLEXINT exited with error code %d.',exitflag);         
    end
    %disp(details.statstring);     
else
    [x,alpha,exitflag]=linprog([1 zeros(1,cols-1)],A,b,Aeq,beq,[zeros(cols_b,1); -inf*ones(K*M,1); zeros(M*Q,1)],[],[],optimset('Display','iter','TolFun',input_parser.Results.TolFun));
    if exitflag~=1
        error('Linprog exited with error code %d.',exitflag); 
    end
end
disp(' ');

F=x(2:2+p-1,1);

end

