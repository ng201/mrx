function [vertices,how]=MRXSEscs(this,T,t)
%MRXSEscs Calculates the extreme vertices of the polyhedral set.
%
% function [vertices,how]=MRXSEscs(this,T,t)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Calculates the extreme vertices of the polyhedral set
%              X = {x: Tx <= t, x >= 0}.               (1)
% Sparse version.
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% T                       - lhs; see (1); sparse matrix
% t                       - rhs; see (1); sparse vector
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% vert                    - vertices
% how                     - result ('ok','failed')
%
%            
% Copyright is with the following author(s):
%
% (C) 2013 Gábor Németh, BME TMIT
%          nemethgab@tmit.bme.hu
% ---------------------------------------------------------------------------
%%
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
%%

switch this.ExtremeSolver

    %==========================================
    % CDD
    %==========================================
    case 5    
        % CDD requires full matrices
        H=struct('A',full(T),'B',full(t));        
        V=cddmex('extreme',H);        
        vertices=sparse(V.V');
        how='ok'; 

    %==========================================
    % lrs
    %==========================================
    case 2
        [vertices,~]=MRXSElrs(T,t);
        how='ok';
        
    %==========================================
    % native matlab function
    %==========================================
    case 1
        
        vertices=MRXSEextrpts(T,t);
        how='ok';
        
    otherwise
  
        %error('MRXSEscf: unknown extreme solver specified!');
        how='failed';
end

