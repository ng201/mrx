function b=contains(pf, theta)
%CONTAINS Report whether the point can be routed in the network.
%
% b=contains(pf, thetab)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% pf                      - pathflow object)
% theta                   - traffic
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% b                       - bool
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

K=pf.K;
P=pf.P;
M=pf.M;
u=pf.u;

p=0;                 % number of path
p_k=zeros(1,K);      % pk(k) is the number of path of k-th sorce-dest
PP=[];
for k=1:K 
    p_k(k)=size(P{k},2);
    p=p+p_k(k);
    PP=[PP P{k}];
end

b=[theta; u];
slacks=M;

% variables u_p \forall p \in P_k, \forall k \in K
A=sparse(M+K,p+slacks);

pi=0;
for k=1:K
    A(k,pi+1:pi+p_k(k))=ones(1,p_k(k));
    pi=pi+p_k(k);
end

pi=0;
for k=1:K
    A(K+1:K+M,pi+1:pi+p_k(k))=P{k};
    pi=pi+p_k(k);
end
A(K+1:K+M,p+1:p+M)=eye(M,M); %add slacks

opts=MRXLoad();
[x,fval,exitflag]=MRXSLscs(opts,zeros(1,p+slacks),[],[],A,b,zeros(p+slacks,1),[],[]);

if exitflag==1
    b=1;    
else
    b=0;
end

end
