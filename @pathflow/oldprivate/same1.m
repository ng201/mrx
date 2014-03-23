function tribool=same1(pf, H, h)
%SAME1 decides whether the polytope given by Hx<h, x>=0 is the same as the thpol of pf.
%
% [T,t]=same1(pf, H, h)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Decides whether the polytope given by Hx<h, x>=0 is the same as the thpol of pf.
% This approach enumerates the extreme points of the polytope Hx<=h, x>=0 and 
% The extreme points of Hx<=h, x>=0 are enumerated and checked one-by-one. 
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% pf                      - pathflow object
% H                       - the lhs of the Hx<=h equation
% h                       - the rhs
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% tribool                 - yes/no/error
%
% (C) 2013 Gábor Németh, BME-TMIT
%          nemethgab@tmit.bme.hu

K=pf.K;

% change almost-zero elements to zero 
H(abs(H)<1e-12)=0;
h(abs(h)<1e-12)=0; 

% test extreme points
P=struct('A',[H;-eye(K)],'B',[h;zeros(K,1)]);

try
    V=cddmex('extreme',P);
    vertices=V.V';
catch
   % CDD sometimes fails to compute extreme points correctly
   % this happens often when slopes of two hyperplanes are too close
   % that's why we use a fixed-point arithmetics
   for ii=1:size(H,1),
       for jj=1:size(H,2),
           H(ii,jj)=fix(H(ii,jj)*(1/1e-5))*1e-5;
       end
       h(ii)=fix(h(ii)*(1/1e-5))*1e-5;
   end  
   [H,h]=reduce(H,h);
   P=struct('A',[H;-eye(K)],'B',[h;zeros(K,1)]);
   %try
       V=cddmex('extreme',P);
       vertices=V.V';
   %end
   
end
        
tribool=1;
for v=1:size(vertices,2)
    b=contains(pf,vertices(:,v));
    tribool=tribool && b;
    if ~tribool
        break;
    end
end

end