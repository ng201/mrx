function this = load(this,varargin)
%LOAD Loads network from the given file.
%
% function this = load(this,varargin)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Load network from the file specified.
%
% Syntax:
%
% g = gio('param1','param2')
%
% load(g,'Directed',true,'Input','fileName','FileType','lgf','CapTag','capacity')
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% this                    - gio object
% varargin                - other parameters, e.g.,
%                           Directed  - whether the network is directed or not
%                           Input     - the name of the file
%                           FileType  - the files type (lgf by default)
%                           CapTag    - the column containing the capacity
%                                       is labeled by CapTag (lgf files, 
%                                       default is cost)
%                           LengthTag - the column containing the capacity
%                                       is labeled by CapTag (lgf files, 
%                                       default is cost)
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% this                    - the modified gio object (note that gio is
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

% parse all the other argument with inputParser
input_parser=inputParser;
input_parser.CaseSensitive = true;

addParameter(input_parser,'Directed',this.directed,@(x) validateattributes(x, {'logical'},{}));
addParameter(input_parser,'Input','',@(x) validateattributes(x, {'char'},{}));
addParameter(input_parser,'Model',this.model,@(x)any(validatestring(x,{'pf','af'})));
addParameter(input_parser,'FileType','.lgf',@(x)any(validatestring(x,{'.i1.j1.k1','lgf'})));
addParameter(input_parser,'CapTag','capacity',@(x) validateattributes(x, {'char'},{}));
addParameter(input_parser,'LengthTag','length',@(x) validateattributes(x, {'char'},{}));

parse(input_parser,varargin{:});

this.directed = input_parser.Results.Directed;
this.model = input_parser.Results.Model;

if ~isempty(find(strcmp('Input',input_parser.UsingDefaults),1))
    err = MException('Input error: Inpur file name cannot be empty');
    throw(err);
end

this.sdpairs=[];
this.P={};

switch input_parser.Results.FileType
    case '.i1.j1.k1'
        this.adj=this.load_i1j1k1(input_parser.Results.Input);
        this.edges=this.generate_el(this.adj);
        this.cap=[];
    case 'lgf'
        [this.adj,this.cap]=this.load_lgf(input_parser.Results.Input,'CapTag',input_parser.Results.CapTag,'LengthTag',input_parser.Results.LengthTag,'Directed',this.directed);
        this.edges=this.generate_el(this.adj);
    otherwise
        
end

end