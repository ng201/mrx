function this = save(this,fname)
%SAVE Saves the gio object to m file.
%
% function this = save(this,fname)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Saves the gio object to an m file.
%
% Syntax:
%
% g = gio('param1','param2')
%
% save(g,'test.m')
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% this                    - gio object
% fname                   - the name of the output file
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% this                    - the gio object (note that gio is
%                           handle)
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

error(nargchk(2,2,nargin));

if isempty(this.adj) 
     error('Network error: please generate or load a new network.');
end
                        
if strcmp(this.model,'pf')                
     if isempty(this.P)
         error('Network error: please define the paths of the path-flow model.');
     end
     
     fid=fopen(fname,'w');
     
     fprintf(fid,'clear P; clear u;\n');
     
     for i=1:length(this.P)
         fprintf(fid,'P{%u} = [',i);
         
         for j=1:length(this.P{i})
             
             fprintf(fid,'%f, %f; ',this.P{i}(j,:));
             
         end
         fprintf(fid,'];\n');
     end
   
     fprintf(fid,'u = [');       
     fprintf(fid,'%f; ',this.edges(:,4));
     fprintf(fid,'];\n');
   
     fclose(fid);
     
else
     error('Not implemented yet.');
end

end