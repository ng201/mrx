function [alpha,F] = oblivious_ratio(pf)
%OBLIVIOUS_RATIO Calculates the oblivious ratio and the oblivious routing function
%
% [alpha,F]=oblivious_ratio(pf, varargin)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Calculates the classic oblivious ratio, i.e., solves the 
%       min alpha : [F_1;...;F_K] = F >= 0, 1^T F_k = 1 k=1,...,K
%                   \forall (i,j) \in E: PFx<=alpha*u
% linear program.
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% pf                      - pathflow object
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% alpha                   - the oblivious ratio
% F                       - the routing function
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

K = pf.K;   % source-dest pairs
M = pf.M;   % number of edges
P = pf.P;   % paths
u = pf.u;   % capacities

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

cols = 1 + p + M*M + K*M;
cols_w= 1 + p;
cols_b= 1 + p + M*M;

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
        A(M+M*p+(e-1)*K+k,cols_b+(e-1)*K+k)=1*u(e);
        A(M+M*p+(e-1)*K+k,pi+1:pi+p_k(k))=P{k}(e,:);
    end
    pi=pi+p_k(k);
end

opts=MRXLoad();
[x,alpha,exitflag]=MRXSLscs(opts,[1 zeros(1,cols-1)],A,b,Aeq,beq,[zeros(1+p+M*M,1); -inf*ones(K*M,1)],[],[]);

if exitflag~=1
     error('The linear program solver exited with error code %d.',exitflag); 
end
disp(' ');

F=x(2:2+p-1,1);

end

