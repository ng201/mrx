function [a,Mf,theta,mflow]=bounding_simplex2(pf, theta0, H, h, varargin)
%BOUNDING_SIMPLEX2 return the equation of the bounding simplex within 
%                  the are given by the cutting planes Hx<=h.
%
% function [M,a,theta]=bounding_simplex2(pf, theta0, H, h [, Options])
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Returns the equation of the bounding simplex in the following form:
%                         ax<=Mf                                       (1)
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% pf                      - pathflow object
% theta0                  - the base poitnt; optional
% H                       - lhs of cutting planes Hx<=h
% h                       - rhs of cutting planes Hx<=h
% Options                 - type 'TryCenter' to place base point in
%                           the middle ('on','off'); default is 'off'
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
% a                       - see equation (1)
% Mf                      - see equation (1)
% theta                   - demand generating exactly M
% mflow                   - the max flows relative to theta0
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

input_parser=inputParser;
input_parser.CaseSensitive = true;

addRequired(input_parser,'pf');

addRequired(input_parser,'theta0',@(x) validateattributes(x, {'numeric'}, {'column','nrows',pf.K}));
addRequired(input_parser,'H',@(x) validateattributes(x, {'numeric'}, {'2d','ncols',pf.K}));
addRequired(input_parser,'h',@(x) validateattributes(x, {'numeric'}, {'column'}));

addParamValue(input_parser,'TryCenter','off',@(x)any(validatestring(x,{'on','off'})));
addParamValue(input_parser,'PositiveOrthant','off',@(x)any(validatestring(x,{'on','off'})));
parse(input_parser,pf,theta0,H,h,varargin{:});

if isempty(theta0)
    theta0=zeros(pf.K,1);
end

if size(H,1)~=length(h)
    error('H must have the same number of rows as h does!');
end

%if size(H,2)~=pf.K
%    error('H should have K columns!');
%end

opts=MRXLoad();

mflow=max_flow(pf,theta0,H,h);

% \sum 1/mf(i) x = Mf + 1/mf * theta0 the equation of the simplex
% Mf == ?

% variables Mf, x[1:K], flow[1:K][1:num_of_flows(i)]
flow_num=0;
for k=1:pf.K
    flow_num=flow_num+size(pf.P{k},2);
end
variables_num=1+pf.K+flow_num;

% objective, we maximixe Mf
f=zeros(1,variables_num);
f(1)=-1;

% EQUALITIES
% setting flows
Aeq=zeros(1+pf.K,variables_num);
pk=1+pf.K;
for k=1:pf.K
    Aeq(k,k+1)=-1;
    
    for l=1:size(pf.P{k},2)
        Aeq(k,pk+l)=1;
    end
    
    pk=pk+size(pf.P{k},2);
end
% setting Mf
Aeq(pf.K+1,1)=-1;
for k=1:pf.K
    Aeq(pf.K+1,1+k)=1/mflow(k);
end
beq=zeros(pf.K+1,1);
beq(pf.K+1)=(1 ./mflow)' * theta0;

% INEQS
% inside TH pol
A=zeros(pf.M,variables_num);
for e=1:pf.M
    pk=1+pf.K;
    for k=1:pf.K
        for l=1:size(pf.P{k},2)
            A(e,pk+l)=pf.P{k}(e,l);
        end
        pk=pk+size(pf.P{k},2);
    end
end
b=pf.u;

% H\theta<=h
AA=[zeros(size(H,1),1) H zeros(size(H,1),flow_num)];
A=[A; AA];
b=[b; h];

% LOWER BOUND
lb=zeros(variables_num,1);
lb(2:1+pf.K)=theta0; % at least theta0

[x,fval,exitflag]=MRXSLscf(opts,f,A,b,Aeq,beq,lb,[],[]);

if exitflag~=1
    error('Linprog exited with error code %d.',exitflag);
end

Mf=-fval+(1 ./mflow)' * theta0;
a=1 ./mflow;
theta=x(2:1+pf.K);


if strcmp(input_parser.Results.TryCenter,'on')
    
    % \alpha_1, ..., \alpha_K, \alpha, f_{1,1}, ..., f_{1,p1}, ...., f_{K,1},
    % ..., f_{K,pK}
    Aeq1=a';
    for i=1:pf.K
        Aeq1(i)=Aeq1(i) * mflow(i);
    end
    
    Aeq1=[Aeq1 0 zeros(1,flow_num); zeros(pf.K,pf.K+1+flow_num)];

    pk=1+pf.K;
    for k=1:pf.K
        
        Aeq1(1+k,k)=mflow(k);
 
        for l=1:size(pf.P{k},2)
            Aeq1(1+k,pk+l)=-1;
        end
    
        pk=pk+size(pf.P{k},2);
    end
    
    beq1=Mf-a' * theta0;
    beq1=[beq1; -theta0];
    
    A1=-eye(pf.K);
    A1=[A1 ones(pf.K,1) zeros(pf.K,flow_num)];
    b1=zeros(pf.K,1);

    A1=[A1; zeros(pf.M,pf.K+1+flow_num)];
    for e=1:pf.M
       pk=1+pf.K;
       for k=1:pf.K
           for l=1:size(pf.P{k},2)
               A1(pf.K+e,pk+l)=pf.P{k}(e,l);
           end
           pk=pk+size(pf.P{k},2);
       end
    end
    b1=[b1; pf.u];
    
    % Hx<=h
    A1=[A1; AA];
    b1=[b1; h];
     
    lb1=zeros(pf.K+1+flow_num,1);
    ub1=[ones(pf.K+1,1); theta0+mflow];

    [x,fval,exitflag]=MRXSLscf(opts,[zeros(1,pf.K) -1 zeros(1,flow_num)],A1,b1,Aeq1,beq1,lb1,ub1,[]);
   
    if exitflag~=1
        error('Linprog exited with error code %d.',exitflag);
    end
    
    theta=x(1:pf.K);
    for i=1:pf.K
        theta(i)=theta(i) * mflow(i)+ theta0(i);
    end

end

if strcmp(input_parser.Results.PositiveOrthant,'on')
    a=[a'; -eye(pf.K)];
    Mf=[Mf; zeros(pf.K,1)];
else
    a=a';
end

end
