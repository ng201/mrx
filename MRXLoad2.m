function mos=MRXLoad2(varargin)
%MRXLOAD Stores and loads an MRX option object (structure).
%
% function mos=MRXLoad(varargin)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Stores and loads an MRX option structure. The value is stored in a
% persistent variable.
%
% Syntax:
%
% mos = MRXLoad('param1','param2')
% mos = MRXLoad
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% varargin                - the arguments in param-value form
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
% mos                     - the mrx options structure
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

persistent pmos;
    
if nargin==1      
    pmos=varargin{1};
end
    
mos=pmos;
    
end