function tribool=same2(pf, H, h)
%SAME1 decides whether the polytope given by Hx<h, x>=0 is the same as the thpol of pf.
%
% [T,t]=same1(pf, H, h, [vertices])
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Decides whether the polytope given by Hx<h, x>=0 is the same as the thpol of pf.
% The extreme points of the dual representation of the THPol are
% enumerated.
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% pf                      - pathflow object
% H                       - the lhs of the Hx<=h equation
% h                       - the rhs
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% tribool                 - yes/no/error
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

narginchk(3,3);

K=pf.K;
P=pf.P;
M=pf.M;
u=pf.u;

p=0;                 % number of path
p_k=zeros(1,K);      % pk(k) is the number of path of k-th sorce-dest
for k=1:K 
    p_k(k)=size(P{k},2);
    p=p+p_k(k);
end

% vars: \beta_k,w_ij
A=sparse(1 + p, K + M);

% first line
A(1,K+1:K+1+M-1)=u';

PP=sparse(p,M);         % paths versus links
PK=sparse(p,K);         % paths versus demand
pi=1;
for k=1:K 
         
    PK(pi:pi+p_k(k)-1,k)=ones(p_k(k),1);
    
    PP(pi:pi+p_k(k)-1,:)=P{k}';
        
    pi=pi+p_k(k);
end

A(2:1+p,1:K)=PK;
A(2:1+p,K+1:K+M)=-PP;

b=sparse(1+p,1);
b(1,1)=1;


opts=MRXLoad();

% test extreme points
P=struct('A',full([A;-eye(K+M)]),'B',full([b;zeros(K+M,1)]));
V=cddmex('extreme',P);
vertices=V.V'
full(vertices(:,1:K))
    
[vertices,how]=MRXSEscs(opts,[A;-eye(K+M)],[b;zeros(K+M,1)]);
%vertices=V.V';
how
full(vertices(:,1:K))

tribool=1;

for v=1:size(vertices,2)
    beta=vertices(:,v)';
    
    %%%%%[x,fval,exitflag]=linprog(-beta(1:K),H,h,[],[],zeros(K,1),[],[],optimset('Display','off'));        
                                      
    [x,fval,exitflag,~]=MRXSLscs(opts,-beta(1:K),H,h,[],[],zeros(K,1),[],[]);
    if exitflag~=1
        error('Linear programm solver exited with error code %d.',exitflag);         
    end
    %disp(details.statstring);
                   
    if -fval>1.0000001
        tribool=0;
        break;
    end
end

end