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