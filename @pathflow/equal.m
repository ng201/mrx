function tribool = equal(pf, H, h)
%EQUAL decides whether the polytope given by Hx<h, x>=0 is the same as the thpol of pf.
%
% tribool = equal(pf, H, h)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Decides whether the polytope given by Hx<h, x>=0 is the same as the thpol of pf.
% This approach uses the Gallo-Ülkücü method. See: G. Gallo, A. Ülkücü, "Bilinear 
% Programming: An Exact Algorithm", Mathematical Programming 12, 1977,
% pp. 173-194
% 
% The orriginal algortihm solves the problem:
%
%             maximize c'*x + x'*Q'*y + d'*y
%                 s.t. A*x <= a,   B'*y <= b
%                        x >= 0,      y >= 0.
% Note that this problem is equivalent to the following peoblem:
%             maximize (c'*x + min b'*u)
%                 s.t. A*x <= a,   B*u >= d + Q*x
%                        X >= 0,     u >= 0.
%
% In the code H and h play the role of A and a, respectively, because we
% can find easily a non-degenerate vertex of Hx<=h, x>=0. This vertex is
% the origin.
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
%
% Copyright is with the following author(s):
%
%
% (C) 2013 Gábor Németh
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

tribool=1;
opts=MRXLoad();

K=pf.K;
P=pf.P;
M=pf.M;
u=pf.u;

R=size(H,2);
% if R=length(h) error!!!!!

p=0;                 % number of path
p_k=zeros(1,K);      % pk(k) is the number of path of k-th sorce-dest
for k=1:K 
    p_k(k)=size(P{k},2);
    p=p+p_k(k);
end

% vars: \beta_k,w_ij
B=sparse(1 + p, K + M);

% first line
B(1,K+1:K+1+M-1)=u';

PP=sparse(p,M);         % paths versus links
PK=sparse(p,K);         % paths versus demand
pi=1;
for k=1:K 
         
    PK(pi:pi+p_k(k)-1,k)=ones(p_k(k),1);
    
    PP(pi:pi+p_k(k)-1,:)=P{k}';
        
    pi=pi+p_k(k);
end

B(2:1+p,1:K)=PK;
B(2:1+p,K+1:K+M)=-PP;

b=sparse(1+p,1);
b(1,1)=1;

%%% H plays the role of A
%%% h plays the role of a

c=sparse(K,1);
d=sparse(K+M,1);
Q=sparse(K+M,K);
for k=1:K
    Q(k,k)=1;
end

% pick a non-degenerate point: origin + the maxflows
V=zeros(K,K+1);
%V=zeros(K,1);
for i=1:K    
    %[theta0,fval,exitflag]=linprog(-1, H(:,i), h, [], [], zeros(1,1));
    theta=zeros(K,1);
    
    theta0b=min(h./H(:,i))-0.000001;
        
    %theta(i,1)=theta0b;
    %V=[V theta];
    %V(:,i+1)=theta;
    V(i,i+1)=theta0b;
end

% step 0.a
zV=0;
for i=1:K 
    x=V(:,i+1);
    zVt=compute_zV(c,Q,d,B',b,x);
    if zVt>zV
        zV=zVt;
        
        % point is outside polytope!
        if zV>1
            tribool=0;
            return;
        end

    end
end


Bhat=[B'; -b'];
dhat=[d; -zV];
Qhat=[Q; c'];

% step 0.b
omega=eye(K);
Delta={};
Delta{1}=omega;

OPTIONS.verbose=0;

while (1)

    % step 1
    if size(Delta,2) == 0 % it is empty
        break;
    end
        
    omega=Delta{1}; % get the first element
    
    % step 2
    % vars theta,u
    theta=zeros(K,1);
    for k=1:K
        v=omega(:,k);
        %[~,fval,exitflag]=linprog([-1 zeros(1,1+p)],[Qhat*v -Bhat],-dhat,[],[],[-inf; zeros(1+p,1)])
        [~,fval,exitflag,~]=MRXSLscs(opts,[-1; zeros(1+p,1)],[Qhat*v -Bhat],-dhat,[],[],[-inf; zeros(1+p,1)],[]);
        if exitflag~=1
            error('111 CPLEXINT exited with error code %d.',exitflag);
        end 
        
        theta(k,1)=-fval;
    end       
    
    % step 3 ???????????
    %[lambda,fval,exitflag]=linprog(-1./theta,H*omega,h,[],[],zeros(K,1))
    [lambda,fval,exitflag,~]=MRXSLscs(opts,(-1./theta),H*omega,h,[],[],zeros(K,1),[]);
    if exitflag~=1
            error('222 CPLEXINT exited with error code %d.',exitflag);
    end    
    
    % step 4
    if -fval<=1.00001
        Delta(1)=[];
        %%%%%warning('Point deleted...');
        continue;
    end
        
    xq=zeros(K,1);
    for k=1:K
        xq=xq+lambda(k)*omega(:,k);
    end
    %%%%%[~,indx]=ismember(bsxfun(@plus,V',-xq')>=0,ones(1,K),'rows');
    [~,indx]=ismember(abs(bsxfun(@plus,V',-xq'))<1e-7,ones(1,K),'rows');
    %%%%%[~,indx]=ismember(V',xq','rows'); % !!!!!!!!!!!!!!!!!!!!!
    
    if find(indx,1)~=0
        Delta(1)=[];
        continue;
    end
    
    V=[V xq];
    zVt=compute_zV(c,Q,d,B',b,xq);
    if zVt>zV
        zV=zVt;
        dhat(K+M+1)=-zV;
        % point is outside polytope!
        if zV>1.0001
            tribool=0;            
            return;
        end
    end
    
    vq=xq/sqrt(sum(xq.^2));
    for k=1:K
        if lambda(k)>0.00001
            
            omega1=omega;
            %omega1(k,:)=vq  %%%%%%%%%%%%%%%%% sor vagy oszlop???
            omega1(:,k)=vq;
            
            Delta{size(Delta,2)+1}=omega1;
            
        end
    end
    Delta(1)=[];

end

end


function zV=compute_zV(c,Q,d,B,b,x)

    %[xx,zV,exitflag]=linprog(b',-B,-d-Q*x,[],[],zeros(size(b,1),1));
    opts=MRXLoad2();    
    [zz,zV,exitflag,~]=MRXSLscs(opts,b,-B,-d-Q*x,[],[],zeros(size(b,1),1),[]);
    if exitflag~=1
        error('zV CPLEXINT exited with error code %d.',exitflag);
    end
    zV=zV+c'*x;
    
end
