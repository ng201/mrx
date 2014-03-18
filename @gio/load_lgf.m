function [adj,cap] = load_lgf(fname,varargin)
%LOAD_LGF loads network topology from file written in lemon graph format.
%
% function [adj,cap] = load_lgf(fname,varargin)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Loads network topology from file written in lemon graph format.
%
% Syntax:
%
% [adj,cap] = load_lgf('textgraph.lgf','CapTag','cap','LengthTag','len','Directed',true);
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% fname                   - file to read
% varargin                - other parameters, e.g.,
%                           Directed  - whether the network is directed or not
%                           CapTag    - the column containing the capacity
%                                       is labeled by CapTag (lgf files, 
%                                       default is cost)
%                           LengthTag - the column containing the capacity
%                                       is labeled by CapTag (lgf files, 
%                                       default is cost)
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% adj                     - the adjacency matrix (with edge lengths)
% cap                     - the matrix of capacities
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

% parse all the other argument with inputParser
input_parser=inputParser;
input_parser.CaseSensitive = true;

addParameter(input_parser,'CapTag','capacity',@(x) validateattributes(x, {'char'},{}));
addParameter(input_parser,'LengthTag','length',@(x) validateattributes(x, {'char'},{}));
addParameter(input_parser,'Directed',false,@(x) validateattributes(x, {'logical'},{}));
parse(input_parser,varargin{:});
    

captag=input_parser.Results.CapTag;
lengthtag=input_parser.Results.LengthTag;

fid=fopen(fname);

tline=fgetl(fid);

while ischar(tline)
    
    if strncmp(tline,'#',1)==1
        tline=fgetl(fid);
        continue;
    end
    
    if strcmp(strtrim(tline),'@arcs')==1 | strcmp(strtrim(tline),'@edges')==1
        break; % trim
    end
    
    tline=fgetl(fid);
end

tline=strtrim(fgetl(fid));
header=textscan(tline,'%s');

capcol=find(ismember(header{1},captag),1)+2
lengthcol=find(ismember(header{1},lengthtag),1)+2

if isempty(capcol) & isempty(lengthcol)
    error('Cap and length not found.');
elseif isempty(capcol)

elseif isempty(lengthcol)
    
end

if capcol < lengthcol
    format = sprintf('%s',['%u' '%u' repmat('%s',1,capcol-1-2) '%f' repmat('%s',lengthcol-capcol-1-2,1) '%f%*[^\n]']);
elseif capcol == lengthcol
    format = sprintf('%s',['%u' '%u' repmat('%s',1,capcol-1-2) '%f%*[^\n]']);
else
    format = sprintf('%s',[repmat('%s',1,lengthcol-1-2) '%f' repmat('%s',capcol-lengthcol-1-2,1) '%f%*[^\n]']);
end

C=textscan(fid,format,'CommentStyle','#');

fclose(fid);

adj=[];
cap=[];
for i=1:length(C{1})
    adj(C{1}(i)+1,C{2}(i)+1)=C{lengthcol}(i);
    cap(C{1}(i)+1,C{2}(i)+1)=C{capcol}(i);
end


% undirected - make it symmetric
if input_parser.Results.Directed == false
    adj=symmetrize(adj);
    cap=symmetrize(cap);
end

end

