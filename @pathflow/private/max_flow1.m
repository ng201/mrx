function mflow=max_flow1(pf, theta0)
%MAX_FLOW1 Calculates max flow
%
% function mf=max_flow1(pf[, theta0])
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Calculates max flows. If theta0 is set, the amount of traffic specified in 
% theta0 has to be routed in the network. The H x <= h inequalities must be true, also. 
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% pf                      - pathflow object
% theta0                  - min. traffic to be routed; optional
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
% mflow                   - the (row vector of) max flows
%
%
% Copyright is with the following author(s):
%
% (C) 2014 G�bor N�meth
%          Inter�University Centre for Telecommunications and Informatics
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

if nargin == 1
    theta0=zeros(pf.K,1);
end

variables_num=1;
for k=1:pf.K            
    variables_num=variables_num+size(pf.P{k},2);
end
    
opts=MRXLoad();
mflow=zeros(pf.K,1);
for k=1:pf.K
        
    % eq
    Aeq=zeros(pf.K,variables_num);
    beq=zeros(pf.K,1);    
    pk=0;
    for kk=1:pf.K
        for l=1:size(pf.P{kk},2)
            Aeq(kk,pk+l)=1;
        end
        if kk~=k            
            beq(kk)=theta0(kk);
        else
            Aeq(kk,variables_num)=-1;
            beq(kk)=0;            
        end
        pk=pk+size(pf.P{kk},2);
    end
   
    % ineq
    A=zeros(pf.M,variables_num);
    for e=1:pf.M
        pk=0;
        for kk=1:pf.K
            for l=1:size(pf.P{kk},2)
                A(e,pk+l)=pf.P{kk}(e,l);
            end
            pk=pk+size(pf.P{kk},2);
        end
    end
    b=pf.u;

    lb=zeros(variables_num,1);
    f=zeros(1,variables_num);
    f(variables_num)=-1;
    
    [x_tmp,fval,exitflag]=MRXSLscf(opts,f,A,b,Aeq,beq,lb,[],[]);

    if exitflag~=1
        error('Linprog exited with error code %d.',exitflag);
    end
    
    mflow(k)=-fval-theta0(k);
    
end

end