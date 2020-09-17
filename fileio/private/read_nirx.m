function varargout = read_nirx(filename, hdr, begsample, endsample, chanindx)

% READ_NIRX reads raw NIRS data from the NIRx wl1/wl2 files
%
% Use as
%   hdr = read_nirx(filename);
%   dat = read_nirx(filename, hdr, begsample, endsample, chanindx);
%   evt = read_nirx(filename, hdr);
%
% The NIRx format consists of a whole bunch of files that go together.
%   filename.avg
%   filename.dat
%   filename.evt
%   filename.hdr
%   filename.inf
%   filename.set
%   filename.tpl
%   filename.wl1
%   filename.wl2
%   filename_config.txt
% This function will read them and represent the combined information in a FieldTYrip
% header structure, event structure or data matrix.
%
% See also FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT, FT_FILETYPE, SNIRF, READ_ARTINIS_OXY3, READ_ARTINIS_OXY4

% Copyright (C) 2020, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

[p, f, x] = fileparts(filename);
filename = fullfile(p, f); % drop the extension

% add the extension for each of the files
avgfile = [filename, '.avg'];
datfile = [filename, '.dat'];
evtfile = [filename, '.evt'];
hdrfile = [filename, '.hdr'];
inffile = [filename, '.inf'];
setfile = [filename, '.set'];
tplfile = [filename, '.tlp'];
wl1file = [filename, '.wl1'];
wl2file = [filename, '.wl2'];
configfile = [filename, '.config.txt'];

needhdr = (nargin==1);
needevt = (nargin==2);
needdat = (nargin==5);

if needhdr
  keyboard
  
elseif needevt
  keyboard
  
elseif needdat
  keyboard
  
end


