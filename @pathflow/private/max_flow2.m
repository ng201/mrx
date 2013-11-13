function mflow=max_flow2(pf, theta0, H, h)
%MAX_FLOW2 Calculates max flow
%
% function mf=max_flow2(pf, theta0, H, h)
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
% H                       - the rhs of the Hx<=h inequality; optional
% h                       - the lhs of the Hx<=h inequality; optional
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
% mflow                   - the (row vector of) max flows
%
% (C) 2012 G�bor N�meth, BME-TMIT
%          nemethgab@tmit.bme.hu

error(nargchk(4,4,nargin));

% if size(H,1)~=length(h)
%     error('H must have the same number of rows as h.');
% end
% 
% if size(H,2)~=pf.K 
%     error('The number of columns in H must be equal with pf.K.');
% end

variables_num=1;
for k=1:pf.K            
    variables_num=variables_num+size(pf.P{k},2);
end

opts=MRXLoad();
mflow=zeros(pf.K,1);
for k=1:pf.K
        
    % eq
    Aeq=zeros(pf.K-1,variables_num);
    beq=zeros(pf.K-1,1);
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
        
    row=size(A,1);
    A=[A; zeros(size(H,1),variables_num)];    
    pk=0;
    for kk=1:pf.K       
        for l=1:size(pf.P{kk},2)
            A(row+1:row+size(H,1),pk+l)=H(:,kk);
        end
        pk=pk+size(pf.P{kk},2);
    end
    b=[b; h];
    
    lb=zeros(variables_num,1);
    f=zeros(variables_num,1);
    f(variables_num)=-1;
        
    [x_tmp,fval,exitflag]=MRXSLscf(opts,f,A,b,Aeq,beq,lb,[],[]);

    if exitflag~=1
        error('Linprog exited with error code %d.',exitflag);
    end
    
    mflow(k)=-fval-theta0(k);
    
end

end