function tribool=same2(pf, H, h, vertices)
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
% vertices                - if it is set, these vertices are used
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% tribool                 - yes/no/error
%
% (C) 2013 Gábor Németh, BME-TMIT
%          nemethgab@tmit.bme.hu

error(nargchk(3,4,nargin));

if nargin == 3
    vertices=[];
end

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

if isempty(vertices)

    % test extreme points
    P=struct('A',full([A;-eye(K+M)]),'B',full([b;zeros(K+M,1)]));
    V=cddmex('extreme',P);
    vertices=V.V';

%     warning off MATLAB:singularMatrix;
%     warning off backtrace;
%     vertices=solvers.extrpts(A,b);        
%     warning(lastwarn);
%     warning on backtrace;
%     warning on MATLAB:singularMatrix;
end

tribool=1;

OPTIONS.verbose=1;
OPTIONS.timelimit=3600;

for v=1:size(vertices,2)
    beta=vertices(:,v);
    
    %%%%[x,fval,exitflag]=linprog(-beta(1:K),H,h,[],[],zeros(K,1),[],[],optimset('Display','off'));
    
    % cplexint

    [x,fval,exitflag,details]=solvers.cplexint.cplexint([],-beta(1:K),H,h,[],[],zeros(K,1),[],[],[],OPTIONS);
    if exitflag~=1
        error('CPLEXINT exited with error code %d.',exitflag);         
    end
    %disp(details.statstring);
                   
    if -fval>1.0000001
        tribool=0;
        break;
    end
end

end