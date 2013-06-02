function cii = read_cifti(filename)

% READ_CIFTI

% 540 bytes with nifti-2 header
% 4 bytes that indicate the presence of a header extension [1 0 0 0]
% 4 bytes with the size of the header extension in big endian?
% 4 bytes with the header extension code NIFTI_ECODE_CIFTI [0 0 0 32]
% variable number of bytes with the xml section, at the end there might be some empty "junk"
% 8 bytes, presumaby with the size and type?
% variable number of bytes with the voxel data

% read the header section
hdr = read_nifti2_hdr(filename);

vox_offset = hdr.vox_offset;
xml_offset = 540+12;
xml_size   = vox_offset-xml_offset-8;

fid = fopen(filename, 'rb', hdr.endian);

fseek(fid, 540, 'bof');
hdrext   = fread(fid, [1 4], 'int8');
if ~isequal(hdrext, [1 0 0 0])
  error('cifti requires a header extension');
end

% this is an alternative to determine the size
prexml = fread(fid, [1 8], 'uint8=>uint8');
if hdr.endian=='l'
  xml_size = typecast(fliplr(prexml(1:4)), 'uint32') - 16;
else
  xml_size = typecast(prexml(1:4), 'uint32') - 8;
end

fseek(fid, xml_offset, 'bof');
xmldata  = fread(fid, [1 xml_size], 'uint8=>char');

if any(xmldata==0)
  warning('removing some junk from xml section');
  xmldata  = xmldata(xmldata>0);
end

% write the xml section to a temporary file
% xmlfile = [tempname, '.xml'];
% fid = fopen(xmlfile, 'w');
% fwrite(fid, xmldata);
% fclose(fid);

% this requires the xmltree object from Guillaume Flandin
% see http://www.artefact.tk/software/matlab/xml/
% it is also included with the gifti toolbox
xml = xmltree(xmldata);


% read the voxel data section
fseek(fid, vox_offset, 'bof');
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

fclose(fid);

cii.hdr = hdr;
cii.xml = xml;
cii.voxdata = squeeze(reshape(voxdata, hdr.dim(2:end)));
cii.xmldata = xmldata;
