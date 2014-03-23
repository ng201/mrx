function [alpha, F, g] = affine_ratio1(af, H, h)
%function [alpha, F, g] = affine_ratio(af, H, h)
%AFFINE_RATIO - Affine oblivious arc-flow routing 
%
% [alpha, F, g] = affine_ratio(af, H, h)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes the affine oblivious routing for the region {t: Ht <= h, t >= 0}
% of the throughput polytope. `Affine' here means that the routing is
% sought for in the form of an affine function: r = Ft + g, t \in H.
% `Oblivious' means that the routing in demand-oblivios: for any throughput
% rquirement in the specified subregion of the throughput polytope, the
% routing si able to accommodate that request in the network with maximally
% `alpha' link-overload . `Arc-flow' routing means that the routing comes 
% in the form of arc-flows (instead of path-flows).
%
% This function solves the giant LP:
%
% min alpha = a
%
% s.t.            NA F_k^k = e_k               \for k
%                 NA F_k^l = 0                 \for k, l \neq k
%                 NA g_k = 0                   \for k
% \for ij \in E:{ (\sum_k g_k^ij)/u_ij + u w^ij + h la^ij <= a
%                 NA^T pi_k^ij <= w^ij        \for k
%                 (e_k)^T pi_k^ij + (H_k)^T la^ij >= 
%                    (\sum_l F_l^ij)^k/u_ij  \for k
%                 w^ij, la^ij >= 0   }
% \for k, ij:{    g_k^ij - v^k,ij - mu^k,ij >= 0
%                 NA^T om_l^k,ij <= v^k,ij     \for k
%                 (e_l)^T om_l^k,ij + (H_l)^T mu^k,ij >= 
%                    -(F_l^ij)^k             \for l
%                 v^k,ij, mu^k,ij >= 0   }
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
%
% NA - node-arc incidence matrix
%
% E - E = { e_k: k \in K }, a matrix whose columns specify the source nodes
% and endpoints of sessions (e_k^i = -1 if i = s_k, e_k^i = 1 if i = d_k,
% e_k^i = 0 othwewise). Again, an arbitrary row is deleted.
% 
% u - capacity vector
%
% H, h - subregion of the throughput polytope for which the oblivious
% routing is computed: {t: Ht <= h, t >= 0) \intersection T(G) (optional,
% a default will be generated if empty)
%
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
%
% alpha - the oblivious ratio
%
% F, g - a cell array containing the affine routing functions, that is, 
% the routing function for the kth session is x_k = F{k} t + g{k}
%
%
% Copyright is with the following author(s):
%
% (C) 2008 Gabor Retvari, BMT-TMIT
%
% (C) 2014 Gábor Németh
%          Inter–University Centre for Telecommunications and Informatics
%          Kassai u. 26., Debrecen, Hungary
%          nemethgab@tmit.bme.hu
%
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

narginchk(1,3);

NA = af.NA;
E = af.E;
u = af.u;

[N, M] = size(NA);

if size(E, 1) ~= N
    error('N, e: size mismatch');
end

if size(u, 1) ~= M
    error('N, u: size mismatch');
end

K = size(E, 2);

% if H is undef
if nargin == 3
    H = -eye(1, K);
    h = 0;
end

if isempty(H) || isempty(h)
    H = -eye(1, K);
    h = 0;
end

Q = size(H, 1);
if size(H, 2) ~= K
    error('E, H: size mismatch');
end

if Q ~= size(h, 1)
    error('H, h: size mismatch');
end

% convert our matrices to sparse
NA = sparse(NA); u = sparse(u); E = sparse(E);
H = sparse(H); h = sparse(h);

cols = 1 + M*K*(1+K) + (K+1)*M*(K*N + M + Q);
Aeq = sparse( N*K*(K+1), cols );
beq = sparse( N*K*(K+1), 1);

disp( sprintf('Problem has\n\t%i columns', cols ));
disp( sprintf('\t%i equalities', N*K*(K+1) ));

% N F_k^l = ...
for k = 1:K
    r = (k-1) * K * N;
    c = (k-1) * K * M;
    e_k = E(:, k);
    for l = 1:K
        rr = r+(l-1)*N;
        cc = c+(l-1)*M;
        Aeq(rr+1:rr+N, cc+1:cc+M) = NA;
        if k == l
            beq(rr+1:rr+N, 1) = e_k;
        end
    end
end

% N g_k = 0
for k = 1:K
    r = K^2 * N + (k-1) * N;
    c = K^2 * M + (k-1) * M;
    for l = 1:K
        Aeq(r+1:r+N, c+1:c+M) = NA;
    end
end

% inequalities
bl_r = 1 + K*M + K;
bl_c = K*N+M+Q;
A = sparse( (K+1)*M*bl_r, cols); 
b = sparse( (K+1)*M*bl_r, 1);

disp( sprintf('\t%i inequalities', (K+1)*M*bl_r ));

% the dual inequalities for \sum_k x_ij_k <= \alpha u_ij 
NT = NA';
for ij = 1:M
    u_ij = u(ij, 1);
    r = (ij - 1)*bl_r;
    c = M*K*(1 + K) + (ij - 1)*bl_c;
    % the first row, difficult one
    r = r + 1;
    % \sum g_k^ij 
    for k = 1:K
        A(r, K^2*M + (k-1) * M + ij) = 1 / u_ij;
    end
    % w^ij u
    A(r, c+K*N+1:c+K*N+M) = u';
    % h la^ij
    A(r, c+K*N+M+1:c+K*N+M+Q) = h';
    % -alpha
    A(r, cols) = -1;

    % N^T pi_k <= w
    for k = 1:K
        A(r+(k-1)*M+1:r+(k-1)*M+M, c+(k-1)*N+1:c+(k-1)*N+N) = NT;
        A(r+(k-1)*M+1:r+(k-1)*M+M, c+K*N+1:c+K*N+M) = ...
            sparse(1:M, 1:M, -1);
    end

    % -(e_k)^T pi_k^ij - (H_k)^T la^ij +
    %           (\sum_l F_l^ij)^k/u_ij <= 0  \for k
    r = r + K*M;
    for k = 1:K
        e_k = E(:, k);
        h_k = H(:, k)';
        A(r+k, c+(k-1)*N+1:c+(k-1)*N+N) = -e_k';
        A(r+k, c+K*N+M+1:c+K*N+M+Q) = -h_k;
        % XXX, check this first, maybe wrong!
        for l = 1:K
           A(r+k, (l-1)*K*M + (k-1)*M + ij) = 1 / u_ij;
        end
    end

    % w^ij, la^ij >= 0
%    r = r + K;
%    A(r+1:r+M+Q, c+K*N+1:c+K*N+M+Q) = ...
%        sparse(1:M+Q, 1:M+Q, -1);
end

% the dual inequalities to x_ij_k >= 0 
for k = 1:K
    for ij = 1:M
        r = M*bl_r + (k-1)*M*bl_r + (ij-1)*bl_r;
        c = M*K*(1+K) + M*bl_c + (k-1)*M*bl_c + (ij-1)*bl_c;
        % the first row, difficult one
        r = r+1;
        % \g_k^ij
        A(r, K^2*M + (k-1) * M + ij) = -1;
        % v^k,ij u
        A(r, c+K*N+1:c+K*N+M) = u';
        % h la^ij
        A(r, c+K*N+M+1:c+K*N+M+Q) = h';
        
        % N^T om_l <= v
        for l = 1:K
            A(r+(l-1)*M+1:r+(l-1)*M+M, c+(l-1)*N+1:c+(l-1)*N+N) = NT;
            A(r+(l-1)*M+1:r+(l-1)*M+M, c+K*N+1:c+K*N+M) = ...
                sparse(1:M, 1:M, -1);
        end
        
        % -(e_l)^T om_l^k,ij - (H_l)^T mu^k,ij +
        %           - F_k^ij)^l <= 0  \for l
        r = r + K*M;
        for l = 1:K
            e_l = E(:, l);
            h_l = H(:, l)';
            A(r+l, c+(l-1)*N+1:c+(l-1)*N+N) = -e_l';
            A(r+l, c+K*N+M+1:c+K*N+M+Q) = -h_l;
            % XXX, check this first, maybe wrong!
            A(r+l, (k-1)*K*M + (l-1)*M + ij) = -1;
        end
        
        % v^k,ij, mu^k,ij >= 0
%        r = r + K;
%        A(r+1:r+M+Q, c+K*N+1:c+K*N+M+Q) = ...
%            sparse(1:M+Q, 1:M+Q, -1);
    end
end

% the objective function
cc = sparse(1, cols);
cc(1, cols) = 1;

% hack up lb
for j = 1:M*K*(1+K)
    lb(j) = -Inf;
end
for ij = 1:M
    c = M*K*(1+K)+(ij-1)*bl_c;
    lb(c+1:c+K*N) = -Inf;
    lb(c+K*N+1:c+K*N+M+Q) = 0;
    
end
for k = 1:K
    for ij = 1:M
        c = M*K*(1+K) + M*bl_c + (k-1)*M*bl_c +(ij-1)*bl_c;
        lb(c+1:c+K*N) = -Inf;
        lb(c+K*N+1:c+K*N+M+Q) = 0;
    end
end
lb(cols) = 0;


opts=MRXLoad();
[x,alpha,exitflag]=MRXSLscs(opts,cc,A,b,Aeq,beq,lb,[],[])

if status ~= 1
    error('The linear program solver exited with error code %d.',exitflag);
end

F = []; g = [];
for k = 1:K
    F_k = zeros(M, K);
    g_k = zeros(M, 1);
    for l = 1:K
        c = (k-1) * K * M + (l-1) * M;
        F_k(1:M, l) = x(c+1:c+M, 1);
    end
    c = K^2 * M + (k-1) * M;
    g_k(1:M, 1) = x(c+1:c+M, 1);
    F{k} = F_k;
    g{k} = g_k;
end

end
