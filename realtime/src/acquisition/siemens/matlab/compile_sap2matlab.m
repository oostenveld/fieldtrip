% SAP2MATLAB parses Siemens ASCII protocol data and generates a 
% corresponding MATLAB data structure.
%
% This function
% is currently used for de-serialising the header information
% from a FieldTrip buffer containing fMRI data from a Siemens
% scanner.

% Copyright (C) 2010, Stefan Klanke,
% 	Modified by Tim van Mourik, 2014
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

warning('Trying to compile MEX file')
oldDir = pwd;
try
  options = '-I../include';
  fileNames = '../src/sap2matlab.c ../src/siemensap.c ';
  libraries = [];
  eval(['mex ' fileNames, libraries, options]);
catch
  rethrow(lasterror)
end

%% Test
load('mrprotString.mat');
S = sap2matlab(apstr);








