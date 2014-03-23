function [H,K] = normalize(H,K,Options)
%NORMALIZE Normalizes a given polytope
%
% [H,h] = normalize(H,h,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% PN = NORMALIZE(P) returns normalized representation PN of a
% polytope P. Polytope PN={x | H_ix <= h_i, i=1,...,nc} has the
% following property: H_i' * H_i=1, for all i.
%
% USAGE:
%   P=normalize(P)
%   P=normalize(P,Options)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P                - Polytope
% Options.rel_tol  - relative tolerance
% Options.abs_tol  - absolute tolerance
%
% Note: If Options is missing or some of the fields are not defined, the 
%       default values will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% [H,h] - polytope in normalized description
% ---------------------------------------------------------------------------

rel_tol = 1e-6;
abs_tol = 1e-7;

[nc,nx]=size(H);
Anorm=sqrt(sum(H .* H,2));

% consistency check
%------------------
if nc<1 | nx<1
    error('NORMALIZE: Polytope has to have nc>0 and nx>0');
end
if any(isinf(Anorm))
   error('NORMALIZE: No Inf terms are allowed in the matrix P.H');
end

% deal with "zero" rows in P.K
%-----------------------------
% use 2-norm combined with absolute and relative tolerance
ii=find(Anorm<=abs_tol | Anorm<=abs(rel_tol*K));

nii=length(ii);
if nii>0
    Anorm(ii)=1;
    H(ii,:)=repmat([1 zeros(1,nx-1)],nii,1);
    
    % decide if some constraint is always true
    jj=(K(ii)>=-abs_tol);
    K(ii(jj))=Inf;
    % or is it always false
    K(ii(~jj))=-Inf;
end

temp=1./Anorm;

H = H .* temp(:,ones(1,nx));   % replaces P.H=P.H .* repmat(temp,1,nx); - call to repmat is rather slow
K = K .* temp;

return
