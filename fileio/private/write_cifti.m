function write_cifti(filename, source, varargin)

% WRITE_CIFTI writes a source structure according to FT_DATATYPE_SOURCE to a cifti file.
%
% Use as
%   write_cifti(filename, source, ...)
% where optional input arguments should come in key-value pairs and may include
%   parameter    = string, fieldname that contains the data
%   parcellation = string, fieldname that descripbes the (optional) parcellation
%
% See also READ_CIFTI

% Copyright (C) 2013, Robert Oostenveld
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
% $Id$

parcellation = ft_getopt(varargin, 'parcellation');
parameter    = ft_getopt(varargin, 'parameter');
precision    = ft_getopt(varargin, 'precision', 'double');

if ~isempty(parcellation)
  assert(ft_datatype(source, 'parcellation') || ft_datatype(source, 'segmentation'), 'the input structure does not define a parcellation');
end

if isfield(source, 'dim')
  % it represents source estimates on regular 3-D grid
  modeltype = 'voxel';
  % the homogenous transformation is required further down
  if ~isfield(source, 'transform')
    source.transform = pos2transform(source.pos);
  end
elseif isfield(source, 'tri')
  % it represents source estimates on a triangulated cortical sheet
  modeltype = 'surface';
else
  % it represents source estimates with an unknown topological arrangement
  modeltype = 'unknown';
end

if isfield(source, 'inside') && islogical(source.inside)
  % convert into an indexed representation
  source.inside = find(source.inside(:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dat = source.(parameter);
if isfield(source, [parameter 'dimord'])
  dimord = source.([parameter 'dimord']);
else
  dimord = source.dimord;
end

dimtok = tokenize(dimord, '_');
if ~strcmp(dimtok{1}, 'pos')
  error('the first dimension should correspond to positions in the brain')
end

switch dimord
  case 'pos_pos'
    if ~isempty(parcellation)
      x = '.ptseries.nii';
    else
      x = '.dtseries.nii';
    end
  case 'pos_time'
    if ~isempty(parcellation)
      x = '.pconn.nii';
    else
      x = '.dconn.nii';
    end
  otherwise
    error('unsupported dimord')
end % switch

[p, f] = fileparts(filename);
filename = fullfile(p, [f x]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct the NIFTI-2 header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if false
  hdr.magic           = fread(fid, [1 8 ], 'int8=>int8'     ); % 4       `n', '+', `2', `\0','\r','\n','\032','\n' or (0x6E,0x2B,0x32,0x00,0x0D,0x0A,0x1A,0x0A)
  hdr.datatype        = fread(fid, [1 1 ], 'int16=>int16'   ); % 12      See file formats
  hdr.dim             = fread(fid, [1 8 ], 'int64=>double'  ); % 16      See file formats
  hdr.vox_offset      = fread(fid, [1 1 ], 'int64=>int64'   ); % 168     Offset of data, minimum=544
  hdr.intent_code     = fread(fid, [1 1 ], 'int32=>int32'   ); % 504     See file formats
  hdr.intent_name     = fread(fid, [1 16], 'int8=>char'     ); % 508     See file formats
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct the XML object describing the geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tree = xmltree;
tree = set(tree, 1, 'name', 'CIFTI');
tree = attributes(tree, 'add', find(tree, 'CIFTI'), 'Version', '1.0');
tree = attributes(tree, 'add', find(tree, 'CIFTI'), 'NumberOfMatrices', '1');
tree = add(tree, find(tree, 'CIFTI'), 'element', 'Matrix');
tree = add(tree, find(tree, 'CIFTI/Matrix'), 'element', 'Volume');

tree = add(tree, find(tree, 'CIFTI/Matrix/Volume'), 'element', 'TransformationMatrixVoxelIndicesIJKtoXYZ');
tree = attributes(tree, 'add', find(tree, 'CIFTI/Matrix/Volume/TransformationMatrixVoxelIndicesIJKtoXYZ'), 'DataSpace', 'NIFTI_XFORM_UNKNOWN');
tree = attributes(tree, 'add', find(tree, 'CIFTI/Matrix/Volume/TransformationMatrixVoxelIndicesIJKtoXYZ'), 'TransformedSpace', 'NIFTI_XFORM_UNKNOWN');
tree = attributes(tree, 'add', find(tree, 'CIFTI/Matrix/Volume/TransformationMatrixVoxelIndicesIJKtoXYZ'), 'UnitsXYZ', 'NIFTI_UNITS_MM');
tree = add(tree, find(tree, 'CIFTI/Matrix/Volume/TransformationMatrixVoxelIndicesIJKtoXYZ'), 'chardata', sprintf('%f ', source.transform));

tree = add(tree, find(tree, 'CIFTI/Matrix'), 'element', 'MatrixIndicesMap');
tree = attributes(tree, 'add', find(tree, 'CIFTI/Matrix/MatrixIndicesMap'), 'IndicesMapToDataType', 'CIFTI_INDEX_TYPE_BRAIN_MODELS');
tree = attributes(tree, 'add', find(tree, 'CIFTI/Matrix/MatrixIndicesMap'), 'AppliesToMatrixDimension', '0');

switch(modeltype)
  case 'voxel'
    
    if ~isempty(parcellation)
      % write one brainmodel per parcel
      pindex = source.([parcellation]);
      plabel = source.([parcellation 'label']);
      for i=1:numel(plabel)
        error('fixme')
      end
      
    elseif isfield(source, 'inside')
      % write one brainmodel, only include the voxels inside the brain
      tree = add(tree, find(tree, 'CIFTI/Matrix/MatrixIndicesMap'), 'element', 'BrainModel');
      tree = attributes(tree, 'add', find(tree, 'CIFTI/Matrix/MatrixIndicesMap/BrainModel'), 'IndexOffset', sprintf('%d ', 0));
      tree = attributes(tree, 'add', find(tree, 'CIFTI/Matrix/MatrixIndicesMap/BrainModel'), 'IndexCount', sprintf('%d ', length(source.inside)));
      tree = attributes(tree, 'add', find(tree, 'CIFTI/Matrix/MatrixIndicesMap/BrainModel'), 'ModelType', 'CIFTI_MODEL_TYPE_VOXELS');
      tree = attributes(tree, 'add', find(tree, 'CIFTI/Matrix/MatrixIndicesMap/BrainModel'), 'BrainStructure', 'CIFTI_STRUCTURE_CORTEX');
      tree = add(tree, find(tree, 'CIFTI/Matrix/MatrixIndicesMap/BrainModel'), 'element', 'VoxelIndicesIJK');
      tree = add(tree, find(tree, 'CIFTI/Matrix/MatrixIndicesMap/BrainModel/VoxelIndicesIJK'), 'chardata', sprintf('%d ', source.inside));
      
    else
      % write one brainmodel for all voxels
      tree = add(tree, find(tree, 'CIFTI/Matrix/MatrixIndicesMap'), 'element', 'BrainModel');
      tree = attributes(tree, 'add', find(tree, 'CIFTI/Matrix/MatrixIndicesMap/BrainModel'), 'IndexOffset', sprintf('%d ', 0));
      tree = attributes(tree, 'add', find(tree, 'CIFTI/Matrix/MatrixIndicesMap/BrainModel'), 'IndexCount', sprintf('%d ', prod(source.dim)));
      tree = attributes(tree, 'add', find(tree, 'CIFTI/Matrix/MatrixIndicesMap/BrainModel'), 'ModelType', 'CIFTI_MODEL_TYPE_VOXELS');
      tree = attributes(tree, 'add', find(tree, 'CIFTI/Matrix/MatrixIndicesMap/BrainModel'), 'BrainStructure', 'CIFTI_STRUCTURE_CORTEX');
      tree = add(tree, find(tree, 'CIFTI/Matrix/MatrixIndicesMap/BrainModel'), 'element', 'VoxelIndicesIJK');
      tree = add(tree, find(tree, 'CIFTI/Matrix/MatrixIndicesMap/BrainModel/VoxelIndicesIJK'), 'chardata', sprintf('%d ', 0:(prod(source.dim-1))));
      
    end % if parcellation, inside or all voxels
    
  case 'surface'
    if ~isempty(parcellation)
      % write one brainmodel per parcel
      pindex = source.([parcellation]);
      plabel = source.([parcellation 'label']);
      for i=1:numel(plabel)
        error('fixme')
      end
      
    elseif isfield(source, 'inside')
      % write one brainmodel, only include the voxels inside the brain
      
    else
      % write one brainmodel for all voxels
      tree = add(tree, find(tree, 'CIFTI/Matrix/MatrixIndicesMap'), 'element', 'BrainModel');
      tree = attributes(tree, 'add', find(tree, 'CIFTI/Matrix/MatrixIndicesMap/BrainModel'), 'IndexOffset', sprintf('%d ', 0));
      tree = attributes(tree, 'add', find(tree, 'CIFTI/Matrix/MatrixIndicesMap/BrainModel'), 'IndexCount', sprintf('%d ', size(source.pos, 1)));
      tree = attributes(tree, 'add', find(tree, 'CIFTI/Matrix/MatrixIndicesMap/BrainModel'), 'ModelType', 'CIFTI_MODEL_TYPE_SURFACE');
      tree = attributes(tree, 'add', find(tree, 'CIFTI/Matrix/MatrixIndicesMap/BrainModel'), 'BrainStructure', 'CIFTI_STRUCTURE_CORTEX');
      tree = attributes(tree, 'add', find(tree, 'CIFTI/Matrix/MatrixIndicesMap/BrainModel'), 'SurfaceNumberOfNodes', sprintf('%d ', size(source.pos, 1)));
      tree = add(tree, find(tree, 'CIFTI/Matrix/MatrixIndicesMap/BrainModel'), 'element', 'NodeIndices');
      tree = add(tree, find(tree, 'CIFTI/Matrix/MatrixIndicesMap/BrainModel/NodeIndices'), 'chardata', sprintf('%d ', 0:(size(source.pos, 1)-1)));
      
    end % if parcellation, inside or all voxels
    
  otherwise
    error('unrecognized description of the geometrical model')
end % case modeltype


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write everything to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 540 bytes with nifti-2 header
% 4 bytes that indicate the presence of a header extension [1 0 0 0]
% 4 bytes with the size of the header extension in big endian?
% 4 bytes with the header extension code NIFTI_ECODE_CIFTI [0 0 0 32]
% variable number of bytes with the xml section, at the end there might be some empty "junk"
% 8 bytes, presumaby with the size and type?
% variable number of bytes with the voxel data

xmlfile = [tempname '.xml'];  % this will contain the cifti XML structure
save(tree, xmlfile);          % write the XMLTREE object to disk

xmlfid = fopen(xmlfile, 'rb');
xmldat = fread(xmlfid, [1, inf], 'char');
fclose(xmlfid);
xmlsize = length(xmldat);
xmlpad  = ceil((xmlsize+8)/16)*16 - (xmlsize+8);

% construct the nifti-2 header
hdr.magic           = [110 43 50 0 13 10 26 10];

% see http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/datatype.html
switch precision
  case 'uint8'
    hdr.datatype = 2;
  case 'int8'
    hdr.datatype = 256;
  case 'uint16'
    hdr.datatype = 512;
  case 'int16'
    hdr.datatype = 4;
  case 'uint32'
    hdr.datatype = 768;
  case 'int32'
    hdr.datatype = 8;
  case 'uint64'
    hdr.datatype = 1280;
  case 'int64'
    hdr.datatype = 1024;
  case 'single'
    hdr.datatype = 16;
  case 'double'
    hdr.datatype = 64;
  otherwise
    error('unsupported precision "%s"', precision);
end

switch precision
  case {'uint8' 'int8'}
    hdr.bitpix = 1*8;
  case {'uint16' 'int16'}
    hdr.bitpix = 2*8;
  case {'uint32' 'int32'}
    hdr.bitpix = 4*8;
  case {'uint64' 'int64'}
    hdr.bitpix = 8*8;
  case 'single'
    hdr.bitpix = 4*8;
  case 'double'
    hdr.bitpix = 8*8;
  otherwise
    error('unsupported precision "%s"', precision);
end

hdr.dim             = [6 0 0 0 0 0 0 0]; % note the 6
hdr.intent_p1       = 0;
hdr.intent_p2       = 0;
hdr.intent_p3       = 0;
hdr.pixdim          = [0 1 1 1 1 1 1 1];
hdr.vox_offset      = 540+8+xmlsize+xmlpad;
hdr.scl_slope       = 0;
hdr.scl_inter       = 0;
hdr.cal_max         = 0;
hdr.cal_min         = 0;
hdr.slice_duration  = 0;
hdr.toffset         = 0;
hdr.slice_start     = 0;
hdr.slice_end       = 0;
hdr.descrip         = char(zeros(1,80));
hdr.aux_file        = char(zeros(1,24));
hdr.qform_code      = 0;
hdr.sform_code      = 0;
hdr.quatern_b       = 0;
hdr.quatern_c       = 0;
hdr.quatern_d       = 0;
hdr.qoffset_x       = 0;
hdr.qoffset_y       = 0;
hdr.qOffset_z       = 0;
hdr.srow_x          = [0 0 0 0];
hdr.srow_y          = [0 0 0 0];
hdr.srow_z          = [0 0 0 0];
hdr.slice_code      = 0;
hdr.xyzt_units      = 0;
hdr.intent_code     = 0;
hdr.intent_name     = char(zeros(1,16));
hdr.dim_info        = 0;
hdr.unused_str      = char(zeros(1,15));

% open the file
fid = fopen(filename, 'wb');

% write the header
write_nifti2_hdr(fid, hdr);

fwrite(fid, [1 0 0 0], 'uint8');
fwrite(fid, 8+xmlsize+xmlpad, 'int32');   % esize
fwrite(fid, 32, 'int32');                 % etype
fwrite(fid, xmldat, 'char');              % write the ascii XML section
fwrite(fid, zeros(1,xmlpad), 'uint8');    % zero-pad to the next 16 byte boundary

fwrite(fid, dat, precision);

fclose(fid);


