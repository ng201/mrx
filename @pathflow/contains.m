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
% (C) 2012 Gábor Németh, BME-TMIT
%          nemethgab@tmit.bme.hu

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
[x,fval,exitflag]=MRXSLscf(opts,zeros(1,p+slacks),[],[],A,b,zeros(p+slacks,1),[],[]);

%INDEQ=[1:1:size(A,1)];
%OPTIONS.verbose=0;
%[x,fval,exitflag,details]=solvers.cplexint.cplexint([],[-ones(p,1); zeros(slacks,1)],A,b,INDEQ,[],zeros(p+slacks,1),[],[],[],OPTIONS);
%[~,~,exitflag,~]=solvers.cplexint.cplexint([],zeros(p+slacks,1),A,b,INDEQ,[],zeros(p+slacks,1),[],[],[],OPTIONS);
%if exitflag~=1
%    error('CPLEXINT exited with error code %d.',exitflag);         
%end
%disp(details.statstring);
%disp(exitflag)
if exitflag==1
    b=1;    
else
    b=0;
end

end
