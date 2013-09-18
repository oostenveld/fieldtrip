function cii = ft_read_cifti(filename, varargin)

% FT_READ_CIFTI reads geometry and functional data or connectivity from a
% cifti file and returns a structure according to FT_DATATYPE_SOURCE.
%
% Use as
%   cii = ft_read_cifti(filename, ...)
% where additional input arguments should be specified as key-value pairs
% and can include
%   representation = ''string', 'tree', 'struct', 'source'
%   readdata       = boolean
%
% See also FT_WRITE_CIFTI, READ_NIFTI2_HDR, WRITE_NIFTI2_HDR

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

representation = ft_getopt(varargin, 'representation', 'source');
geometry       = ft_getopt(varargin, 'geometry', {'midthickness', 'pial', 'white', 'inflated', 'very_inflated', 'sphere'});
readdata       = ft_getopt(varargin, 'readdata', false);

% convert 'yes'/'no' into boolean
readdata = istrue(readdata);

if ~iscell(geometry)
  geometry = {geometry};
end

% read the header section
hdr = read_nifti2_hdr(filename);

% xml_offset = 540+12;
% xml_size   = hdr.vox_offset-xml_offset-8;

fid = fopen(filename, 'rb', hdr.endian);

% determine the file size, this is used to catch endian errors
fseek(fid, 0, 'eof');
filesize = ftell(fid);
fseek(fid, 0, 'bof');

fseek(fid, 540, 'bof');
hdrext = fread(fid, [1 4], 'int8');
if hdrext(1)~=1
  error('cifti requires a header extension');
end

% determine the size of the header extension
esize = fread(fid, 1, 'int32=>int32');
etype = fread(fid, 1, 'int32=>int32');

hdrsize = 540;
voxsize = filesize-hdr.vox_offset;
if esize>(filesize-hdrsize-voxsize)
  warning('the endianness of the header extension is inconsistent with the nifti-2 header');
  esize = swapbytes(esize);
  etype = swapbytes(etype);
end

if etype~=32 && etype~=swapbytes(int32(32)) % FIXME there is an endian problem
  error('the header extension type is not cifti');
end

% read the extension content, subtract the 8 bytes from esize and etype
xmldata = fread(fid, [1 esize-8], 'uint8=>char');

% the size of the extension must be an integer multiple of 16 bytes according to http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/extension.html
% consequently there might be some zero-padding at the end of the XML section
if any(xmldata==0)
  xmldata = xmldata(xmldata>0);
end

% write the xml section to a temporary file
xmlfile = 'test.xml';
tmp = fopen(xmlfile, 'w');
fwrite(tmp, xmldata);
fclose(tmp);

% this requires the xmltree object from Guillaume Flandin
% see http://www.artefact.tk/software/matlab/xml/
% it is also included with the gifti toolbox

switch representation
  case 'char'
    cii.xmldata = xmldata;
  case 'tree'
    cii.tree = xmltree(xmldata);
  case 'struct'
    tmp = xmltree(xmldata);
    cii = tree2struct(tmp);
  case 'source'
    tmp = xmltree(xmldata);
    cii = tree2struct(tmp); % the rest of the conversion is done further down
  otherwise
    error('unsupported representation')
end

% include the nifti-2 header in the output structure
cii.hdr = hdr;

if readdata
  % read the voxel data section
  fseek(fid, hdr.vox_offset, 'bof');
  switch hdr.datatype
    case   2, [voxdata nitemsread] = fread(fid, inf, 'uchar');
    case   4, [voxdata nitemsread] = fread(fid, inf, 'short');
    case   8, [voxdata nitemsread] = fread(fid, inf, 'int');
    case  16, [voxdata nitemsread] = fread(fid, inf, 'float');
    case  64, [voxdata nitemsread] = fread(fid, inf, 'double');
    case 512, [voxdata nitemsread] = fread(fid, inf, 'ushort');
    case 768, [voxdata nitemsread] = fread(fid, inf, 'uint');
    otherwise, error('unsupported datatype');
  end
  cii.data = squeeze(reshape(voxdata, hdr.dim(2:end)));
end
fclose(fid);

% try to get the geometrical information from a corresponding gifti files
% the following assumec the convention of the Human Connectome Project
[p, f, x] = fileparts(filename);
t = tokenize(f, '.');

if length(t)==4
  subject  = t{1};
  dataname = t{2};
  geomodel = t{3};
  content  = t{4};
elseif length(t)==5
  subject  = t{1};
  dataname = [t{2} '.' t{3}];
  geomodel = t{4};
  content  = t{5};
else
  warning('cannot decipher the file name, not reading geometry');
  geometry = {};
end

for i=1:length(geometry)
  Lfilename = fullfile(p, [subject '.L.' geometry{i} '.' geomodel '.surf.gii']);
  Rfilename = fullfile(p, [subject '.R.' geometry{i} '.' geomodel '.surf.gii']);
  if exist(Lfilename, 'file') && exist(Rfilename, 'file')
    mesh = ft_read_headshape({Lfilename, Rfilename});
    cii.pos = mesh.pnt;
    cii.tri = mesh.tri;
    break % only read a single mesh
  end
end

if strcmp(representation, 'source')
  cii = struct2source(cii);
end

if readdata
  % rename the data field
  cii.(fixname(dataname)) = cii.data;
  cii = rmfield(cii, 'data');
  
  % rename the datalabel field
  if isfield(cii, 'datalabel')
    cii.(fixname([dataname 'label'])) = cii.datalabel;
    cii = rmfield(cii, 'datalabel');
  end
end

return % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the cifti XML section can be represented in various manners
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cifti = tree2struct(tree)

numericAttributeTypes = {'NumberOfMatrices', 'AppliesToMatrixDimension', 'IndexOffset', 'IndexCount', 'SurfaceNumberOfNodes', 'VolumeDimensions'};

Cifti            = struct(); % the parent of the XML tree, it only contains version info
MatrixIndicesMap = struct(); % this is the interesting content

attr = attributes(tree, 'get', 1);
for j=1:length(attr)
  if any(strcmp(attr{j}.key, numericAttributeTypes))
    Cifti.(attr{j}.key) = str2num(attr{j}.val);
  else
    Cifti.(attr{j}.key) = attr{j}.val;
  end
end

uid_Volume = find(tree,'/CIFTI/Matrix/Volume');
if ~isempty(uid_Volume)
  volume = branch(tree, uid_Volume);
  attr = attributes(volume, 'get', 1); % there is only one attribute here
  if any(strcmp(attr.key, numericAttributeTypes))
    Volume.(attr.key) = str2num(attr.val);
  else
    Volume.(attr.key) = attr.val;
  end
  uid_Transform = find(volume,'/Volume/TransformationMatrixVoxelIndicesIJKtoXYZ');
  transform = branch(volume, uid_Transform);
  attr = attributes(transform, 'get', 1);
  for j=1:length(attr)
    if any(strcmp(attr{j}.key, numericAttributeTypes))
      Volume.(attr{j}.key) = str2num(attr{j}.val);
    else
      Volume.(attr{j}.key) = attr{j}.val;
    end
  end
  Volume.Transform = str2num(get(transform, 2, 'value'));
  Volume.Transform = reshape(Volume.Transform, [4 4])'; % it needs to be transposed
else
  Volume = [];
end

uid_MatrixIndicesMap = find(tree,'/CIFTI/Matrix/MatrixIndicesMap');
for i=1:length(uid_MatrixIndicesMap)
  map = branch(tree, uid_MatrixIndicesMap(i));
  
  % get the attributes of each map
  attr = attributes(map, 'get', 1);
  for j=1:length(attr)
    if any(strcmp(attr{j}.key, numericAttributeTypes))
      MatrixIndicesMap(i).(attr{j}.key) = str2num(attr{j}.val);
    else
      MatrixIndicesMap(i).(attr{j}.key) = attr{j}.val;
    end
  end
  
  uid_NamedMap = find(map, '/MatrixIndicesMap/NamedMap');
  for j=1:length(uid_NamedMap)
    namedmap = branch(map, uid_NamedMap(j));
    MatrixIndicesMap(i).MapName = get(namedmap, children(namedmap, find(namedmap, '/NamedMap/MapName')), 'value');
    uid_LabelTable = find(namedmap, '/NamedMap/LabelTable');
    for k=1:length(uid_LabelTable);
      labeltable = branch(namedmap, uid_LabelTable(k));
      uid_Label = find(labeltable, '/LabelTable/Label');
      for l=1:length(uid_Label)
        % there are also potentially intersting atributes here, but I don't know what to do with them
        MatrixIndicesMap(i).NamedMap.LabelTable.Label{l} = get(labeltable, children(labeltable, uid_Label(l)), 'value');
        attr = attributes(branch(labeltable, uid_Label(l)), 'get', 1);
        for m=1:numel(attr)
          switch attr{m}.key
            case 'Key'
              MatrixIndicesMap(i).NamedMap.LabelTable.Key(l)   = str2double(attr{m}.val);
            case 'Red'
              MatrixIndicesMap(i).NamedMap.LabelTable.Red(l)   = str2double(attr{m}.val);
            case 'Green'
              MatrixIndicesMap(i).NamedMap.LabelTable.Green(l) = str2double(attr{m}.val);
            case 'Blue'
              MatrixIndicesMap(i).NamedMap.LabelTable.Blue(l)  = str2double(attr{m}.val);
            case 'Alpha'
              MatrixIndicesMap(i).NamedMap.LabelTable.Alpha(l) = str2double(attr{m}.val);
          end
        end
      end
    end
  end
  
  uid_BrainModel = find(map, '/MatrixIndicesMap/BrainModel');
  for j=1:length(uid_BrainModel)
    brainmodel = branch(map, uid_BrainModel(j));
    
    % get the attributes of each model
    attr = attributes(brainmodel, 'get', 1);
    for k=1:length(attr)
      if any(strcmp(attr{k}.key, numericAttributeTypes))
        MatrixIndicesMap(i).BrainModel(j).(attr{k}.key) = str2num(attr{k}.val);
      else
        MatrixIndicesMap(i).BrainModel(j).(attr{k}.key) = attr{k}.val;
      end
    end % for
    
    switch MatrixIndicesMap(i).BrainModel(j).ModelType
      case 'CIFTI_MODEL_TYPE_SURFACE'
        uid = find(brainmodel, '/BrainModel/NodeIndices');
        if ~isempty(uid)
          MatrixIndicesMap(i).BrainModel(j).NodeIndices = str2num(get(brainmodel, children(brainmodel, uid), 'value'));
        elseif MatrixIndicesMap(i).BrainModel(j).IndexCount==MatrixIndicesMap(i).BrainModel(j).SurfaceNumberOfNodes
          % assume that there is a one-to-one mapping
          MatrixIndicesMap(i).BrainModel(j).NodeIndices = (1:MatrixIndicesMap(i).BrainModel(j).IndexCount)-1;
        end
        
      case 'CIFTI_MODEL_TYPE_VOXELS'
        MatrixIndicesMap(i).BrainModel(j).VoxelIndicesIJK = str2num(get(brainmodel, children(brainmodel, find(brainmodel, '/BrainModel/VoxelIndicesIJK')), 'value'));
        
      otherwise
        error('unsupported ModelType');
    end % switch
  end % for each BrainModel
end % for each MatrixIndicesMap

% add it to the main structure
Cifti.MatrixIndicesMap = MatrixIndicesMap;
Cifti.Volume           = Volume;
return % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the cifti XML section can be represented in various manners
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function source = struct2source(Cifti)

MatrixIndicesMap = Cifti.MatrixIndicesMap;
dimord = {};

for i=1:length(MatrixIndicesMap)
  switch MatrixIndicesMap(i).IndicesMapToDataType
    case 'CIFTI_INDEX_TYPE_BRAIN_MODELS'
      dimord(MatrixIndicesMap(i).AppliesToMatrixDimension+1) = {'pos'};
      
      TotalNumberOfNodes  = 0;
      TotalNumberOfValues = 0;
      
      % concatenate all graynode positions
      for j=1:length(MatrixIndicesMap(i).BrainModel)
        posbeg = TotalNumberOfValues + 1;
        posend = TotalNumberOfValues + MatrixIndicesMap(i).BrainModel(j).IndexCount;
        
        % in the cifti-1 specification the IndexOffset is not what you
        % would expect it to be, hence we have to use the TotalNumberOfNodes
        % dataIndex(posbeg:posend) = MatrixIndicesMap(i).BrainModel(j).NodeIndices + MatrixIndicesMap(i).BrainModel(j).IndexOffset + 1;
        
        dataIndex(posbeg:posend) = MatrixIndicesMap(i).BrainModel(j).NodeIndices + TotalNumberOfNodes + 1;
        TotalNumberOfNodes  = TotalNumberOfNodes  + MatrixIndicesMap(i).BrainModel(j).SurfaceNumberOfNodes;
        TotalNumberOfValues = TotalNumberOfValues + MatrixIndicesMap(i).BrainModel(j).IndexCount;
        
        ModelType(posbeg:posend) = j; % indexed representation, see ft_datatype_parcellation
        ModelTypelabel{j} = MatrixIndicesMap(i).BrainModel(j).ModelType;
        
        BrainStructure(posbeg:posend) = j; % indexed representation, see ft_datatype_parcellation
        BrainStructurelabel{j} = MatrixIndicesMap(i).BrainModel(j).BrainStructure;
        
      end % for
      
      % it is nicer to have them as column vector
      ModelType      = ModelType';
      BrainStructure = BrainStructure';
      
    case 'CIFTI_INDEX_TYPE_SCALARS'
      dimord{MatrixIndicesMap(i).AppliesToMatrixDimension+1} = []; % scalars are not explicitly represented
      
    case 'CIFTI_INDEX_TYPE_LABELS'
      key = MatrixIndicesMap(i).NamedMap.LabelTable.Key;
      lab = MatrixIndicesMap(i).NamedMap.LabelTable.Label;
      sel = key>0;
      Cifti.datalabel(key(sel)) = lab(sel);
      
    case 'CIFTI_INDEX_TYPE_ TIME_POINTS'
      dimord(MatrixIndicesMap(i).AppliesToMatrixDimension+1) = {'time'};
      keyboard
      
    case 'CIFTI_INDEX_TYPE_FIBERS'
      error('not yet implemented');
      
    case 'CIFTI_INDEX_TYPE_PARCELS'
      error('not yet implemented');
      
    otherwise
      error('unsupported IndicesMapToDataType');
  end % switch
end

dimord = dimord(~isempty(dimord));
source.dimord = sprintf('%s_', dimord{:});
source.dimord(end) = [];

source.pos = Cifti.pos;
source.tri = Cifti.tri;
Nvertices  = size(source.pos,1);
Nvalues    = length(dataIndex);

if isfield(Cifti, 'data')
  % make the data consistent with the graynode positions
  switch source.dimord
    case 'pos'
      tmp = nan(Nvertices,1);
      tmp(dataIndex,:) = Cifti.data;
      
      % case 'pos_time'
      % case 'pos_pos'
      % case 'pos_pos_time'
    otherwise
      error('not yet implemented');
  end % switch
  
  source.data = tmp;
  if isfield(Cifti, 'datalabel')
    source.datalabel = Cifti.datalabel;
  end
end % if data

tmp = nan(Nvertices,1);
tmp(dataIndex,:) = ModelType;
source.ModelType = tmp;
source.ModelTypelabel = ModelTypelabel;

tmp = nan(Nvertices,1);
tmp(dataIndex,:) = BrainStructure;
source.BrainStructure = tmp;
source.BrainStructurelabel = BrainStructurelabel;

return % function

