function cii = ft_read_cifti(filename, varargin)

% FT_READ_CIFTI reads geometry and functional data or connectivity from a
% cifti file and returns a structure according to FT_DATATYPE_SOURCE.
%
% Use as
%   cii = ft_read_cifti(filename, ...)
% where additional input arguments should be specified as key-value pairs
% and can include
%   representation = string, can be 'tree', 'struct' or 'source'
%   readdata       = boolean, can be false or true
%
% See also FT_WRITE_CIFTI, READ_NIFTI2_HDR, WRITE_NIFTI2_HDR

% Copyright (C) 2013-2014, Robert Oostenveld
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
readdata       = ft_getopt(varargin, 'readdata', []); % default depends on file size

% convert 'yes'/'no' into boolean
readdata = istrue(readdata);

% read the header section
hdr = read_nifti2_hdr(filename);

% xml_offset = 540+12;
% xml_size   = hdr.vox_offset-xml_offset-8;

fid = fopen(filename, 'rb', hdr.endian);

% determine the file size, this is used to catch endian errors
fseek(fid, 0, 'eof');
filesize = ftell(fid);
fseek(fid, 0, 'bof');

% set the default for readdata
if isempty(readdata)
  if filesize>1e9
    warning('filesize>1GB, not reading data by default. Please specify ''readdata'' option.');
    readdata = false;
  else
    readdata = true;
  end
end

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
    cii.tree = xmltree(xmldata); % using xmltree from artefact.dk
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
  fprintf('reading grayordinate data...\n');
  switch hdr.datatype
    case   2, [voxdata, nitemsread] = fread(fid, inf, 'uchar');   assert(nitemsread>0);
    case   4, [voxdata, nitemsread] = fread(fid, inf, 'short');   assert(nitemsread>0);
    case   8, [voxdata, nitemsread] = fread(fid, inf, 'int');     assert(nitemsread>0);
    case  16, [voxdata, nitemsread] = fread(fid, inf, 'float');   assert(nitemsread>0);
    case  64, [voxdata, nitemsread] = fread(fid, inf, 'double');  assert(nitemsread>0);
    case 512, [voxdata, nitemsread] = fread(fid, inf, 'ushort');  assert(nitemsread>0);
    case 768, [voxdata, nitemsread] = fread(fid, inf, 'uint');    assert(nitemsread>0);
    otherwise, error('unsupported datatype');
  end
  fprintf('finished reading grayordinate data\n');
  cii.data = squeeze(reshape(voxdata, hdr.dim(2:end)));
end
fclose(fid);


if strcmp(representation, 'source')
  cii = struct2source(cii);
end

% try to get the geometrical information from a corresponding gifti files
% the following assumec the convention of the Human Connectome Project
[p, f, x] = fileparts(filename);
t = tokenize(f, '.');

subject  = 'unknown';
dataname = 'unknown';
geomodel = '';

if length(t)==2
  subject  = t{1};
  dataname = t{2};
elseif length(t)==3
  subject  = t{1};
  dataname = t{2};
  content  = t{3};
elseif length(t)==4
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
  error('cannot parse file name');
end

% the surface anatomy is represented in an external file
% which can be difficult to match with the data

Lfilelist = {
  [subject '.L' '.midthickness'  '.' geomodel '.surf.gii']
  [subject '.L' '.pial'          '.' geomodel '.surf.gii']
  [subject '.L' '.white'         '.' geomodel '.surf.gii']
  [subject '.L' '.inflated'      '.' geomodel '.surf.gii']
  [subject '.L' '.very_inflated' '.' geomodel '.surf.gii']
  [subject '.L' '.sphere'        '.' geomodel '.surf.gii']
  [subject '.L' '.'              '.' geomodel '.surf.gii']
  [subject '.L' '.midthickness'               '.surf.gii']
  [subject '.L' '.pial'                       '.surf.gii']
  [subject '.L' '.white'                      '.surf.gii']
  [subject '.L' '.inflated'                   '.surf.gii']
  [subject '.L' '.very_inflated'              '.surf.gii']
  [subject '.L' '.sphere'                     '.surf.gii']
  [subject '.L'                               '.surf.gii']
  [subject '.CIFTI_STRUCTURE_CORTEX_LEFT'     '.surf.gii']
  };

Rfilelist = {
  [subject '.R' '.midthickness'  '.' geomodel  '.surf.gii']
  [subject '.R' '.pial'          '.' geomodel  '.surf.gii']
  [subject '.R' '.white'         '.' geomodel  '.surf.gii']
  [subject '.R' '.inflated'      '.' geomodel  '.surf.gii']
  [subject '.R' '.very_inflated' '.' geomodel  '.surf.gii']
  [subject '.R' '.sphere'        '.' geomodel  '.surf.gii']
  [subject '.R' '.'              '.' geomodel  '.surf.gii']
  [subject '.R' '.midthickness'                '.surf.gii']
  [subject '.R' '.pial'                        '.surf.gii']
  [subject '.R' '.white'                       '.surf.gii']
  [subject '.R' '.inflated'                    '.surf.gii']
  [subject '.R' '.very_inflated'               '.surf.gii']
  [subject '.R' '.sphere'                      '.surf.gii']
  [subject '.R'                                '.surf.gii']
  [subject '.CIFTI_STRUCTURE_CORTEX_RIGHT'     '.surf.gii']
  };

Bfilelist = {
  [subject '.midthickness'  '.' geomodel '.surf.gii']
  [subject '.pial'          '.' geomodel '.surf.gii']
  [subject '.white'         '.' geomodel '.surf.gii']
  [subject '.inflated'      '.' geomodel '.surf.gii']
  [subject '.very_inflated' '.' geomodel '.surf.gii']
  [subject '.sphere'        '.' geomodel '.surf.gii']
  [subject                  '.' geomodel '.surf.gii']
  [subject '.midthickness'               '.surf.gii']
  [subject '.pial'                       '.surf.gii']
  [subject '.white'                      '.surf.gii']
  [subject '.inflated'                   '.surf.gii']
  [subject '.very_inflated'              '.surf.gii']
  [subject '.sphere'                     '.surf.gii']
  [subject                               '.surf.gii']
  [subject '.CIFTI_STRUCTURE_CORTEX'     '.surf.gii']
  };


if all(ismember({'CIFTI_STRUCTURE_CORTEX_LEFT', 'CIFTI_STRUCTURE_CORTEX_RIGHT'}, cii.BrainStructurelabel))
  for i=1:length(Lfilelist)
    Lfilename = fullfile(p, Lfilelist{i});
    Rfilename = fullfile(p, Rfilelist{i});
    
    if exist(Lfilename, 'file') && exist(Rfilename, 'file')
      warning('reading left hemisphere geometry from %s',  Lfilename);
      meshL = ft_read_headshape(Lfilename);
      warning('reading right hemisphere geometry from %s',  Rfilename);
      meshR = ft_read_headshape(Rfilename);
      
      indexL = find(cii.BrainStructure==find(strcmp(cii.BrainStructurelabel, 'CIFTI_STRUCTURE_CORTEX_LEFT')));
      indexR = find(cii.BrainStructure==find(strcmp(cii.BrainStructurelabel, 'CIFTI_STRUCTURE_CORTEX_RIGHT')));
      
      cii.pos(indexL,:) = meshL.pnt;
      cii.pos(indexR,:) = meshR.pnt;
      
      cii.tri = [
        indexL(meshL.tri)
        indexR(meshR.tri)
        ];
      
      break % only read a single pair of meshes
    end
  end
elseif ismember({'CIFTI_STRUCTURE_CORTEX'}, cii.BrainStructurelabel)
  for i=1:length(Bfilelist)
    Bfilename = fullfile(p, Bfilelist{i});
    
    if exist(Bfilename, 'file')
      warning('reading surface geometry from %s',  Bfilename);
      meshB   = ft_read_headshape(Bfilename);
      indexB  = find(cii.BrainStructure==find(strcmp(cii.BrainStructurelabel, 'CIFTI_STRUCTURE_CORTEX')));
      cii.pos(indexB,:) = meshB.pnt;
      cii.tri = indexB(meshB.tri);
    end
    
    break % only read a single mesh
  end
end



if readdata
  if isfield(cii, 'data')
    % rename the data field
    cii.(fixname(dataname)) = cii.data;
    cii = rmfield(cii, 'data');
  end
  
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

numericAttributeTypes = {'NumberOfMatrices', 'AppliesToMatrixDimension', 'IndexOffset', 'IndexCount', 'SurfaceNumberOfNodes', 'VolumeDimensions', 'SurfaceNumberOfVertices', 'SeriesStart', 'SeriesStep', 'NumberOfSeriesPoints', 'SeriesExponent'};

Cifti            = struct(); % the parent of the XML tree, it only contains version info
MatrixIndicesMap = struct(); % this is the interesting content

attr = attributes(tree, 'get', 1);
if ~iscell(attr), attr = {attr}; end % treat one attribute just like multiple attributes
for j=1:length(attr)
  if any(strcmp(attr{j}.key, numericAttributeTypes))
    Cifti.(attr{j}.key) = str2num(attr{j}.val);
  else
    Cifti.(attr{j}.key) = attr{j}.val;
  end
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
  
  uid_Volume = find(tree,'/CIFTI/Matrix/MatrixIndicesMap/Volume');
  % the following will fail if there are multiple volumes
  if ~isempty(uid_Volume)
    volume = branch(tree, uid_Volume);
    attr = attributes(volume, 'get', 1); % there should only be one attribute here
    if ~iscell(attr), attr = {attr}; end % treat one attribute just like multiple attributes
    for j=1:length(attr)
      if any(strcmp(attr{j}.key, numericAttributeTypes))
        Volume.(attr{j}.key) = str2num(attr{j}.val);
      else
        Volume.(attr{j}.key) = attr{j}.val;
      end
    end
    uid_Transform = find(volume,'/Volume/TransformationMatrixVoxelIndicesIJKtoXYZ');
    if ~isempty(uid_Transform)
      transform = branch(volume, uid_Transform);
      attr = attributes(transform, 'get', 1);
      if isstruct(attr), attr = {attr}; end % treat one attribute just like multiple attributes
      for j=1:length(attr)
        if any(strcmp(attr{j}.key, numericAttributeTypes))
          Volume.(attr{j}.key) = str2num(attr{j}.val);
        else
          Volume.(attr{j}.key) = attr{j}.val;
        end
      end
      Volume.Transform = str2num(get(transform, 2, 'value'));
      Volume.Transform = reshape(Volume.Transform, [4 4])'; % it needs to be transposed
    end
  else
    Volume = [];
  end
  
  uid_NamedMap = find(map, '/MatrixIndicesMap/NamedMap');
  for j=1:length(uid_NamedMap)
    namedmap = branch(map, uid_NamedMap(j));
    MatrixIndicesMap(i).NamedMap(j).MapName = get(namedmap, children(namedmap, find(namedmap, '/NamedMap/MapName')), 'value');
    uid_LabelTable = find(namedmap, '/NamedMap/LabelTable');
    for k=1:length(uid_LabelTable);
      labeltable = branch(namedmap, uid_LabelTable(k));
      uid_Label = find(labeltable, '/LabelTable/Label');
      for l=1:length(uid_Label)
        % there are also potentially intersting atributes here, but I don't know what to do with them
        MatrixIndicesMap(i).NamedMap(j).LabelTable.Label{l} = get(labeltable, children(labeltable, uid_Label(l)), 'value');
        attr = attributes(branch(labeltable, uid_Label(l)), 'get', 1);
        for m=1:numel(attr)
          switch attr{m}.key
            case 'Key'
              MatrixIndicesMap(i).NamedMap(j).LabelTable.Key(l)   = str2double(attr{m}.val);
            case 'Red'
              MatrixIndicesMap(i).NamedMap(j).LabelTable.Red(l)   = str2double(attr{m}.val);
            case 'Green'
              MatrixIndicesMap(i).NamedMap(j).LabelTable.Green(l) = str2double(attr{m}.val);
            case 'Blue'
              MatrixIndicesMap(i).NamedMap(j).LabelTable.Blue(l)  = str2double(attr{m}.val);
            case 'Alpha'
              MatrixIndicesMap(i).NamedMap(j).LabelTable.Alpha(l) = str2double(attr{m}.val);
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
        switch Cifti.Version
          case {'1' '1.0'}
            uid = find(brainmodel, '/BrainModel/NodeIndices');
            try, MatrixIndicesMap(i).BrainModel(j).NodeIndices = str2num(get(brainmodel, children(brainmodel, uid), 'value')); end
          case {'2' '2.0'}
            uid = find(brainmodel, '/BrainModel/VertexIndices');
            try, MatrixIndicesMap(i).BrainModel(j).VertexIndices = str2num(get(brainmodel, children(brainmodel, uid), 'value')); end
          otherwise
            error('unsupported version');
        end % switch version
        
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
dimord = cell(size(Cifti.MatrixIndicesMap));

for i=1:length(MatrixIndicesMap)
  switch MatrixIndicesMap(i).IndicesMapToDataType
    case 'CIFTI_INDEX_TYPE_BRAIN_MODELS'
      dimord(MatrixIndicesMap(i).AppliesToMatrixDimension+1) = {'pos'};
      
      IndexOffset         = [MatrixIndicesMap(i).BrainModel(:).IndexOffset];
      IndexCount          = [MatrixIndicesMap(i).BrainModel(:).IndexCount];
      ModelType           = nan(length(MatrixIndicesMap(i).BrainModel),1);
      ModelTypelabel      = cell(length(MatrixIndicesMap(i).BrainModel),1);
      BrainStructure      = nan(length(MatrixIndicesMap(i).BrainModel),1);
      BrainStructurelabel = cell(length(MatrixIndicesMap(i).BrainModel),1);
      
      tmp = cumsum([0 IndexCount]);
      if ~isequal(IndexOffset, tmp(1:end-1))
        % this happens in some of the example cifti1 files
        % and might be a bug in the actual format of the data in those files
        warning('inconsistency between IndexOffset and IndexCount');
      end
      
      % count the number of greynodes in all combined brain models
      geomCount = 0;
      Cifti.pos = nan(0,3);
      
      % concatenate all greynode positions
      for j=1:length(MatrixIndicesMap(i).BrainModel)
        
        switch MatrixIndicesMap(i).BrainModel(j).ModelType
          case 'CIFTI_MODEL_TYPE_SURFACE'
            switch Cifti.Version
              case {'1' '1.0'}
                posbeg = geomCount + 1;
                posend = geomCount + MatrixIndicesMap(i).BrainModel(j).SurfaceNumberOfNodes;
                geomCount = geomCount + MatrixIndicesMap(i).BrainModel(j).SurfaceNumberOfNodes; % increment with the number of vertices in the (external) surface
                Cifti.pos(posbeg:posend,:) = nan;
                Cifti.dataIndex{j}         = IndexOffset(j) + (1:IndexCount(j));  % these are indices in the data
                Cifti.greynodeIndex{j}     = posbeg:posend;                       % these are indices in the greynodes (vertices or subcortical voxels)
                if isfield(MatrixIndicesMap(i).BrainModel(j), 'NodeIndices')
                  % data is only present on a subset of vertices
                  Cifti.greynodeIndex{j} = Cifti.greynodeIndex{j}(MatrixIndicesMap(i).BrainModel(j).NodeIndices+1);
                end
                
              case {'2' '2.0'}
                posbeg = geomCount + 1;
                posend = geomCount + MatrixIndicesMap(i).BrainModel(j).SurfaceNumberOfVertices;
                geomCount = geomCount + MatrixIndicesMap(i).BrainModel(j).SurfaceNumberOfVertices; % increment with the number of vertices in the (external) surface
                Cifti.pos(posbeg:posend,:) = nan;
                Cifti.dataIndex{j}         = IndexOffset(j) + (1:IndexCount(j));  % these are indices in the data
                Cifti.greynodeIndex{j}     = posbeg:posend;                       % these are indices in the greynodes (vertices or subcortical voxels)
                if isfield(MatrixIndicesMap(i).BrainModel(j), 'VertexIndices')
                  % data is only present on a subset of vertices
                  Cifti.greynodeIndex{j} = Cifti.greynodeIndex{j}(MatrixIndicesMap(i).BrainModel(j).VertexIndices+1);
                end
                
              otherwise
                error('unsupported version');
            end % switch version
            
            
          case 'CIFTI_MODEL_TYPE_VOXELS'
            posbeg = geomCount + 1;
            posend = geomCount + IndexCount(j);
            geomCount = geomCount + IndexCount(j); % increment with the number of vertices in the subcortical structure
            
            Cifti.pos(posbeg:posend,:) = reshape(MatrixIndicesMap(i).BrainModel(j).VoxelIndicesIJK, 3, IndexCount(j))';
            Cifti.dataIndex{j}         = IndexOffset(j) + (1:IndexCount(j));  % these are indices in the data
            Cifti.greynodeIndex{j}     = posbeg:posend;                       % these are indices in the greynodes (vertices or subcortical voxels)
            
          otherwise
            error('unexpected ModelType');
        end % switch
        
        % perform a sanity check on the data and greynode indices
        assert(numel(Cifti.dataIndex{j})==numel(Cifti.greynodeIndex{j}));
        
        ModelType(posbeg:posend) = j; % indexed representation, see ft_datatype_parcellation
        ModelTypelabel{j} = MatrixIndicesMap(i).BrainModel(j).ModelType;
        
        BrainStructure(posbeg:posend) = j; % indexed representation, see ft_datatype_parcellation
        BrainStructurelabel{j} = MatrixIndicesMap(i).BrainModel(j).BrainStructure;
        
      end % for all BrainModels
      
    case 'CIFTI_INDEX_TYPE_SCALARS'
      dimord{MatrixIndicesMap(i).AppliesToMatrixDimension+1} = []; % scalars are not explicitly represented
      if isfield(MatrixIndicesMap(i), 'NamedMap')
        for j=1:length(MatrixIndicesMap(i).NamedMap)
          Cifti.mapname{j} = fixname(MatrixIndicesMap(i).NamedMap(j).MapName);
        end
      end
      
    case 'CIFTI_INDEX_TYPE_LABELS'
      dimord{MatrixIndicesMap(i).AppliesToMatrixDimension+1} = []; % labels are not explicitly represented
      for j=1:length(MatrixIndicesMap(i).NamedMap)
        key = MatrixIndicesMap(i).NamedMap(j).LabelTable.Key;
        lab = MatrixIndicesMap(i).NamedMap(j).LabelTable.Label;
        sel = key>0;
        Cifti.labeltable{j}(key(sel)) = lab(sel);
        Cifti.mapname{j} = fixname(MatrixIndicesMap(i).NamedMap(j).MapName);
      end
      
    case 'CIFTI_INDEX_TYPE_SERIES'
      % this only applies to cifti version 2
      switch MatrixIndicesMap(i).SeriesUnit
        case 'SECOND'
          dimord(MatrixIndicesMap(i).AppliesToMatrixDimension+1) = {'time'};
          Cifti.time = (((1:MatrixIndicesMap(i).NumberOfSeriesPoints)-1) * MatrixIndicesMap(i).SeriesStep + MatrixIndicesMap(i).SeriesStart) * 10^MatrixIndicesMap(i).SeriesExponent;
        case 'HZ'
          dimord(MatrixIndicesMap(i).AppliesToMatrixDimension+1) = {'freq'};
          Cifti.freq = (((1:MatrixIndicesMap(i).NumberOfSeriesPoints)-1) * MatrixIndicesMap(i).SeriesStep + MatrixIndicesMap(i).SeriesStart) * 10^MatrixIndicesMap(i).SeriesExponent;
        otherwise
          error('unsupported SeriesUnit');
      end % switch
      
      
    case 'CIFTI_INDEX_TYPE_TIME_POINTS'
      % this only applies to cifti version 1
      dimord(MatrixIndicesMap(i).AppliesToMatrixDimension+1) = {'time'};
      switch MatrixIndicesMap(i).TimeStepUnits
        case 'NIFTI_UNITS_SEC'
          Cifti.fsample = 1/str2double(MatrixIndicesMap(i).TimeStep);
        otherwise
          % other units should be trivial to implement
          error('unsupported TimeStepUnits');
      end
      
      % case 'CIFTI_INDEX_TYPE_PARCELS'
      %   error('not yet implemented');
      
      % case 'CIFTI_INDEX_TYPE_FIBERS'
      %   error('not yet implemented');
      
    otherwise
      error('unsupported IndicesMapToDataType');
  end % switch
end

dimord = dimord(~cellfun(@isempty, dimord));
source.dimord = sprintf('%s_', dimord{:});
source.dimord(end) = [];

source.pos = Cifti.pos;
Ngreynodes = size(source.pos,1);

if issubfield(Cifti, 'Volume.Transform')
  % this only applies to the voxel coordinates, not to surface vertices which are NaN
  source.pos = ft_warp_apply(Cifti.Volume.Transform, source.pos);
end

if isfield(Cifti, 'data')
  % make the data consistent with the graynode positions
  dataIndex     = [Cifti.dataIndex{:}];
  greynodeIndex = [Cifti.greynodeIndex{:}];
  
  switch source.dimord
    case 'pos'
      [m, n] = size(Cifti.data);
      if m>n
        dat = nan(Ngreynodes,n);
        dat(greynodeIndex(dataIndex),1) = Cifti.data;
      else
        dat = nan(Ngreynodes,m);
        dat(greynodeIndex(dataIndex),:) = transpose(Cifti.data);
      end
    case 'pos_pos'
      dat = nan(Ngreynodes,Ngreynodes);
      dat(greynodeIndex(dataIndex),greynodeIndex(dataIndex)) = Cifti.data;
    case 'pos_time'
      Ntime = size(Cifti.data,2);
      dat = nan(Ngreynodes,Ntime);
      dat(greynodeIndex(dataIndex),:) = Cifti.data;
    case 'time_pos'
      Ntime = size(Cifti.data,1);
      dat = nan(Ngreynodes,Ntime);
      dat(greynodeIndex(dataIndex),:) = transpose(Cifti.data);
      source.dimord = 'pos_time';
    otherwise
      error('unsupported dimord');
  end % switch
  
  if isfield(Cifti, 'mapname') && length(Cifti.mapname)>1
    % use distict names if there are multiple scalars or labels
    for i=1:length(Cifti.mapname)
      fieldname = fixname(Cifti.mapname{i});
      source.(fieldname) = dat(:,i);
      if isfield(Cifti, 'labeltable')
        source.([fieldname 'label']) = Cifti.labeltable{i};
      end
    end
  else
    % the name of the data will be based on the filename
    source.data = dat;
  end
end % if data

source = copyfields(Cifti, source, {'time', 'freq'});

if ~isempty(Cifti.Volume)
  source.dim        = Cifti.Volume.VolumeDimensions;
  source.transform  = Cifti.Volume.Transform;
end

source.BrainStructure      = BrainStructure;
source.BrainStructurelabel = BrainStructurelabel;

return % function
