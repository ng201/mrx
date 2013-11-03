function out=sI(sclass,stype)
%MPT_SI returns information about a given solver
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%
% internal function
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
%%
% Copyright is with the following author(s):
%
% (C) 2013 Gábor Németh, BME-TMIT
%          nemethgab@tmit.bme.hu

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
%%

error(nargchk(2,2,nargin));

if ischar(stype),
    
    if strcmpi(sclass, 'lp')
        convfunc = @char2int_lpsolver;
    elseif strcmpi(sclass, 'extreme')
        convfunc = @char2int_extreme;
    else
        error(sprintf('unknwon solver class ''%s'' !', sclass));
    end
else    
    if strcmpi(sclass, 'lp')
        convfunc = @int2char_lpsolver;
    elseif strcmpi(sclass, 'extreme')
        convfunc = @int2char_extreme;
    else
        error(sprintf('unknwon solver class ''%s'' !', sclass));
    end
end
    
out = feval(convfunc,stype);

return

end


%------------------------------------------------------------------------
function out = int2char_lpsolver(solver);

switch solver
    
    case 0 
                                        out = 'auto';
    case 1 
                                        out = 'linprog';
    case 2 
                                        out = 'glpkmex';
    case 3 
                                        out = 'cplexint';
    case 4 
                                        out = 'cplex';
    case 5 
                                        out = 'cdd criss-cross';
    case 6 
                                        out = 'cdd dual-simplex';
    case 7 
                                        out = 'MOSEK';
    case 901
                                        out = 'MPT2';
    case 902
                                        out = 'MPT3';
    otherwise
                                        out = 'unknown'; 
end

end


%------------------------------------------------------------------------
function out = int2char_extreme(solver);

switch solver
    
    case 0 
                                out = 'auto';
    case 1
                                out = 'extrpts';
    case 2 
                                out = 'lrs';
    case 5 
                                out = 'CDD';
    otherwise
                                out = 'unknown'; 
end

end


%------------------------------------------------------------------------
function out = char2int_lpsolver(solver);
% converts a string description to numerical values

switch solver
             
    case 'linprog'
                                                out=1;
    case {'glpk','glpkmex','glpkcc'}
                                                out=2;
    case 'cplexint'
                                                out=3;
    case {'cplex','Cplex'}  
                                                out=4;
    case {'cdd.cris-cross','cdd.criss-cros', ...
          'cdd.cris-cros','cdd.criss-cross', ...
          'criss-cross','cris-cross', ...
          'criss-cros','cris-cros', ...
          'cdd criss-cross'}
                                                out=5;
    case {'cdd.dual-simplex','cdd dual-simplex', ...
          'cdd dual simplex'}
                                                out=6;
    case 'mosek'
                                                out=7;
               
    case 'mpt2'
                                                out=901;
    case 'mpt'
                                                out=902;  
    otherwise
                                                out=-1;
end

end


%------------------------------------------------------------------------
function out = char2int_extreme(solver);
% converts a string description to numerical values

switch solver
    case {'cdd','CDD'}
                               out=5;
    case 'extrpts' 
                               out=1;
    case 'lrs'
                               out=2;                                   
    otherwise
                               out=0;
end

end

