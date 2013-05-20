function [hdr] = read_nifti2_hdr(filename)

% READ_NIFTI2_HDR

fid = fopen(filename, 'rb', 'l');
hdr.sizeof_hdr = fread(fid, [1 1 ], 'int32=>int32'   ); % 1

% this is borrowed from the FreeSurfer nifti-1 implementation
if hdr.sizeof_hdr~=348 && hdr.sizeof_hdr~=540
  fclose(fid);
  % Now try opening as big endian
  fid = fopen(filename, 'r', 'b');
  hdr.sizeof_hdr = fread(fid, [1 1 ], 'int32=>int32'   ); % 1
  if hdr.sizeof_hdr~=348 && hdr.sizeof_hdr~=540
    fclose(fid);
    error('cannot open %s as nifti file, hdr size = %d, should be 348 or 540\n', filename, hdr.sizeof_hdr);
  end
  hdr.endian = 'l';
else
  hdr.endian = 'b';
end

if hdr.sizeof_hdr==384
  % the remainder of the code is for nifti-2 files
  error('%s seems to be a nifti-1 file', filename)
end

hdr.magic           = fread(fid, [1 8 ], 'int8=>int8'     ); % 4       `n', '+', `2', `\0','\r','\n','\032','\n' or (0x6E,0x2B,0x32,0x00,0x0D,0x0A,0x1A,0x0A)
hdr.datatype        = fread(fid, [1 1 ], 'int16=>int16'   ); % 12      See file formats
hdr.bitpix          = fread(fid, [1 1 ], 'int16=>int16'   ); % 14      See file formats
hdr.dim             = fread(fid, [1 8 ], 'int64=>double'  ); % 16      See file formats
hdr.intent_p1       = fread(fid, [1 1 ], 'double=>double' ); % 80      0
hdr.intent_p2       = fread(fid, [1 1 ], 'double=>double' ); % 88      0
hdr.intent_p3       = fread(fid, [1 1 ], 'double=>double' ); % 96      0
hdr.pixdim          = fread(fid, [1 8 ], 'double=>double' ); % 104     0,1,1,1,1,1,1,1
hdr.vox_offset      = fread(fid, [1 1 ], 'int64=>int64'   ); % 168     Offset of data, minimum=544
hdr.scl_slope       = fread(fid, [1 1 ], 'double=>double' ); % 176     1
hdr.scl_inter       = fread(fid, [1 1 ], 'double=>double' ); % 184     0
hdr.cal_max         = fread(fid, [1 1 ], 'double=>double' ); % 192     0
hdr.cal_min         = fread(fid, [1 1 ], 'double=>double' ); % 200     0
hdr.slice_duration  = fread(fid, [1 1 ], 'double=>double' ); % 208     0
hdr.toffset         = fread(fid, [1 1 ], 'double=>double' ); % 216     0
hdr.slice_start     = fread(fid, [1 1 ], 'int64=>int64'   ); % 224     0
hdr.slice_end       = fread(fid, [1 1 ], 'int64=>int64'   ); % 232     0
hdr.descrip         = fread(fid, [1 80], 'int8=>char'     ); % 240     All zeros
hdr.aux_file        = fread(fid, [1 24], 'int8=>char'     ); % 320     All zeros
hdr.qform_code      = fread(fid, [1 1 ], 'int32=>int32'   ); % 344     NIFTI_XFORM_UNKNOWN (0)
hdr.sform_code      = fread(fid, [1 1 ], 'int32=>int32'   ); % 348     NIFTI_XFORM_UNKNOWN (0)
hdr.quatern_b       = fread(fid, [1 1 ], 'double=>double' ); % 352     0
hdr.quatern_c       = fread(fid, [1 1 ], 'double=>double' ); % 360     0
hdr.quatern_d       = fread(fid, [1 1 ], 'double=>double' ); % 368     0
hdr.qoffset_x       = fread(fid, [1 1 ], 'double=>double' ); % 376     0
hdr.qoffset_y       = fread(fid, [1 1 ], 'double=>double' ); % 384     0
hdr.qOffset_z       = fread(fid, [1 1 ], 'double=>double' ); % 392     0
hdr.srow_x          = fread(fid, [1 4 ], 'double=>double' ); % 400     0,0,0,0
hdr.srow_y          = fread(fid, [1 4 ], 'double=>double' ); % 432     0,0,0,0
hdr.srow_z          = fread(fid, [1 4 ], 'double=>double' ); % 464     0,0,0,0
hdr.slice_code      = fread(fid, [1 1 ], 'int32=>int32'   ); % 496     0
hdr.xyzt_units      = fread(fid, [1 1 ], 'int32=>int32'   ); % 500     0xC (seconds, millimeters)
hdr.intent_code     = fread(fid, [1 1 ], 'int32=>int32'   ); % 504     See file formats
hdr.intent_name     = fread(fid, [1 16], 'int8=>char'     ); % 508     See file formats
hdr.dim_info        = fread(fid, [1 1 ], 'int8=>int8'     ); % 524     0
hdr.unused_str      = fread(fid, [1 15], 'int8=>char'     ); % 525     All zeros
% disp(ftell(fid));                                          % 540     End of the header

fclose(fid);