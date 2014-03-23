function [H,K,keptrows]=reduce(H,K,Options,lambda)
%REDUCE Reduces the polytope by removing redundant inequalities
%
% [H,h,keptrows] = reduce)(H,h,Options,lambda)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% [P,keptrows]=REDUCE(P,OPTIONS) computes the minimal representation of the 
% polytope P by eliminating all redundant constraints. The reduction of P is 
% implemented by solving one LP for each facet. 
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------  
% P                     - Polytope
% Options.lpsolver      - LP solver to be used (see help mpt_solveLP)
% Options.abs_tol       - absolute tolerance
% Options.calculatebox  - Use a bounding box to discard certain hyperplanes a 
%                         priori (1: use BoundingBox)
% lambda                - Optional: Set of lagrange multipliers from chebyball; 
%                         If available, this speeds up the reduction procedure.
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% P         - polytope in reduced description
% keptrows  - indices of hyperplanes which are non-redundand
%
% see also POLYTOPE, NORMALIZE
%

% ---------------------------------------------------------------------------
% Copyright is with the following author(s):
%
% (C) 2004 Raphael Suard, Automatic Control Laboratory, ETH Zurich,
%          suardr@control.ee.ethz.ch
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%          baotic@control.ee.ethz.ch

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
% --------------------------------------------------------------------------

if nargin<3,
    Options=[];
end
if nargin<4,
    lambda=[];
end

rel_tol = 1e-6;
abs_tol = 1e-7;

lpsolver=3;

keptrows = [];

[H,K]=normalize(H,K);
Horig = H;
Korig = K;

% If any boundary is -Inf polytope P is empty
%--------------------------------------------
if any(K==-Inf)
    n=size(H,2);
    H = [1 zeros(1,n-1)];
    K = -Inf;
    return;
end

% Remove rows with Inf boundaries
%--------------------------------
ii=(K==Inf);
H(ii,:)=[];
K(ii)=[];

[nc,nx] = size(H);

if nc==0
    H=[1 zeros(1,nx-1)];
    K=Inf;
    error('REDUCE: Polytope P = R^nx');     % If P is full space R^nx it is an error
end

xCheb = [];
RCheb = [];
if isempty(xCheb),
    H = H;
    K = K;
    Opt.lpsolver = 3;
    [xCheb,RCheb,lambda]=chebyball(H,K,Opt);
end

if RCheb<abs_tol
    H = 1;
    K = -Inf;
    return;
end


%--------------------------------------------------------------------------
%Use a bounding box to discard certain hyperplanes a priori

% Initialize Values
In=eye(nx);
l=zeros(nx,1);               % Lower bounds
u=zeros(nx,1);               % Upper bounds

% Build Objective Matrices

% Determine external box ; minimize x_i
for i=1:nx,
    f=In(i,:);
    
    x=mpt_solveLPs(f,H,K,[],[],xCheb,lpsolver);
    
    l(i)=x(i);
end

% Determine external box ; maximize x_i
for i=1:nx,
    f=-In(i,:);
    
    x=mpt_solveLPs(f,H,K,[],[],xCheb,lpsolver);
    
    u(i)=x(i);
end

% Eliminating HP out the BB
cand = ~((Horig>0).*Horig*(u-l) - (Korig-Horig*l) < -max(1e-4,abs_tol));

keptrows = find(cand); %keep only the candidates
H=Horig(keptrows,:);
K=Korig(keptrows);
nc = numel(K);
P.bbox = [l u];
%--------------------------------------------------------------------------


% Eliminate redundances otherwise same hP two times in P => keptrows
% error

[a,b] = unique([H K],'rows');
b = sort(b);        % we must sort, otherwise we don't respect order of hyperplanes
keptrows = keptrows(b);
H = Horig(keptrows, :);
K = Korig(keptrows, :);
nc = length(K);
cand=ones(nc,1);


%--------------------------------------------------------------------------
% hit and run
x0 = xCheb;
n = nx;
r = nc/2;
tol10 = 10*abs_tol;
keep = zeros(1,2*r);

for i = 1:r
    d = randn(n,1);
    d = d/norm(d);
    den = H*d;
    if any(den==0),
        continue
    end
    t = (K-H*x0)./(H*d);
    pos = t>0;
    negpos = ~pos;
    bound_pos = find(pos);
    bound_neg = find(negpos);
    tt = t(bound_pos);
    [t_closest_pos,j1] = min(tt);
    if isempty(j1),
        continue
    end
    if length(find(abs(t-t_closest_pos)<tol10))>1,
        continue
    end
    tt = t(bound_neg);
    [t_closest_neg,j2] = max(tt);
    if isempty(j2),
        continue
    end
    if length(find(abs(t-t_closest_neg)<tol10))>1,
        continue
    end
    if ~isempty(j1),
        keep(2*i-1)= bound_pos(j1);
    else
        t_closest_pos = 0;
    end
    if ~isempty(j2),
        keep(2*i)  = bound_neg(j2);
    else
        t_closest_neg = 0;
    end
    x1 = x0+d*t_closest_pos;
    x2 = x0+d*t_closest_neg;
    x0 = x1+(x2-x1)/2;
end

keep = unique(keep);

% if "ismembc" is not found on your system, uncomment line below
cand=(~ismembc(1:length(cand), keep))';
%cand=(~ismember(1:length(cand), keep))';

%--------------------------------------------------------------------------

if nc==1
    return;
else

    % Compute minimal representation
    %-------------------------------

%         if (lpsolver==3 | lpsolver==5)
%             % switched off in this version, feature subject to testing
%             %%%%%%%%%%%%%%%%%%%%%%%%%%
%             % use fast call to cddmex
%             %=========================
%     
%             Q.A=P.H;
%             Q.B=P.K;
%             [Pr,killed]=cddmex('reduce_h',Q);          % reduce H representation
%             P.H=Pr.A;
%             P.K=Pr.B;
%             keptrows = setdiff(1:length(Q.B),killed);
%     
%             return;
%         end

    % initialize set of non-redundant hyperplanes
    nonredundant=logical(ones(nc,1));  
    removerow=[];

    k=1;
    while k<=length(cand)
        f_cand = find(cand);
        if any(k==f_cand)
            f1 = H(k,:);
            K(k) = K(k)+0.1;
            
            [xopt,fval,lambda,exitflag,status]=mpt_solveLPs(-f1,H(nonredundant,:),K(nonredundant),[],[],xCheb,lpsolver);
            
            K(k) = K(k)-0.1;
            obj=f1*xopt-K(k);
            nonredundant(k)=0;
            if obj>abs_tol,
                % non-redundant
                nonredundant(k)=1; 
            elseif strcmp(status,'unbounded'),
                % non-redundant
                nonredundant(k)=1;
            else
                removerow=[removerow k];
            end
        end
        k=k+1;
    end

    keptrows(removerow)=[];
    keptrows=sort(keptrows);
    H(~nonredundant,:)=[];
    K(~nonredundant)=[];
    H = H;
    K = K;

    if isempty(keptrows)
        % no constraints => polytope is R^n
        error('REDUCE: Polytope P = R^nx');
    end
end
