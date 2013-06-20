function cii = read_cifti(filename)

% READ_CIFTI reads geometry and functional data or connectivity from a cifti file
%
% Use as
%   cii = read_cifti(filename)
%
% See also WRITE_CIFTI

% 540 bytes with nifti-2 header
% 4 bytes that indicate the presence of a header extension [1 0 0 0]
% 4 bytes with the size of the header extension in big endian?
% 4 bytes with the header extension code NIFTI_ECODE_CIFTI [0 0 0 32]
% variable number of bytes with the xml section, at the end there might be some empty "junk"
% 8 bytes, presumaby with the size and type?
% variable number of bytes with the voxel data

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

hdrsize = 540+4;
voxsize = filesize-hdr.vox_offset;
if esize>(filesize-hdrsize-voxsize)
  warning('the endianness of the header extension is inconsistent with the nifti-2 header');
  esize = swapbytes(esize);
  etype = swapbytes(etype);
end

if etype~=32
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
% xmlfile = [tempname, '.xml'];
% fid = fopen(xmlfile, 'w');
% fwrite(fid, xmldata);
% fclose(fid);

% this requires the xmltree object from Guillaume Flandin
% see http://www.artefact.tk/software/matlab/xml/
% it is also included with the gifti toolbox
xml = xmltree(xmldata);

if false
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
end

fclose(fid);

cii.hdr = hdr;
cii.xml = xml;
% cii.voxdata = squeeze(reshape(voxdata, hdr.dim(2:end)));
cii.xmldata = xmldata;
