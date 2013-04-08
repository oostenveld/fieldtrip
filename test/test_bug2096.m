function test_bug2096

% TEST test_bug2096
% TEST ft_sourcewrite

load(dccnfilename('/home/common/matlab/fieldtrip/data/test/bug2096/CP10168_4DEXP_3-Restin_BNN_V1_MEG_icaimagcoh_freq3.mat'));

cfg = [];
cfg.filename = [tempname '.cii'];
ft_sourcewrite(cfg, source);

% test reading these files
% http://brainvis.wustl.edu/cifti/DenseConnectome.dconn.nii
% http://brainvis.wustl.edu/cifti/DenseTimeSeries.dtseries.nii
% http://brainvis.wustl.edu/cifti/ParcellatedTimeSeries.ptseries.nii

p = '/home/common/matlab/fieldtrip/data/test/bug2096';
% cii1 = ft_read_mri(fullfile(p, 'DenseConnectome.dconn.nii',          'fileformat', 'cifti');
% cii2 = ft_read_mri(fullfile(p, 'DenseTimeSeries.dtseries.nii',       'fileformat', 'cifti');
% cii3 = ft_read_mri(fullfile(p, 'ParcellatedTimeSeries.ptseries.nii', 'fileformat', 'cifti');
cii1 = read_cifti(fullfile(p, 'DenseConnectome.dconn.nii');
cii2 = read_cifti(fullfile(p, 'DenseTimeSeries.dtseries.nii');
cii3 = read_cifti(fullfile(p, 'ParcellatedTimeSeries.ptseries.nii');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cii = read_cifti(filename)
% read the header section
hdr = read_nifti2_hdr(filename);
xml_offset = 540;
dat_offset = hdr.vox_offset;
% read the xml section
fid = fopen(filename, 'rb');
fseek(fid, xml_offset, bof);
buf = fread(fid, [1 dat_offset-xml_offset, 'uint8=>uint8');
% read the voxel data section
switch hdr.datatype
 case   2, [hdr.vol nitemsread] = fread(fp,inf,'uchar');
 case   4, [hdr.vol nitemsread] = fread(fp,inf,'short');
 case   8, [hdr.vol nitemsread] = fread(fp,inf,'int');
 case  16, [hdr.vol nitemsread] = fread(fp,inf,'float');
 case  64, [hdr.vol nitemsread] = fread(fp,inf,'double');
 case 512, [hdr.vol nitemsread] = fread(fp,inf,'ushort');
 case 768, [hdr.vol nitemsread] = fread(fp,inf,'uint');
 otherwise, error('unsupported datatype');
end
fclose(fid);
% conveft the XML section, this requires an external toolbox
ft_hastoolbox('xml4mat', 1);
xml = xml4mat(buf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hdr = read_nifti2_hdr(filename)
fid = fopen(filename, 'rb', 'ieee-le');
hdr.magic           = fread(fid, [1 8 ], 'int8=>int8'     ); % 4       `n', '+', `2', `\0','\r','\n','\032','\n' or (0x6E,0x2B,0x32,0x00,0x0D,0x0A,0x1A,0x0A)
hdr.datatype        = fread(fid, [1 1 ], 'int16=>int16'   ); % 12      See file formats
hdr.bitpix          = fread(fid, [1 1 ], 'int16=>int16'   ); % 14      See file formats
hdr.dim             = fread(fid, [1 8 ], 'int64=>int64'   ); % 16      See file formats
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
hdr.descrip         = fread(fid, [1 80], 'int8=>int8'     ); % 240     All zeros
hdr.aux_file        = fread(fid, [1 24], 'int8=>int8'     ); % 320     All zeros
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
hdr.unused_str      = fread(fid, [1 15], 'int8=>int8'     ); % 525     All zeros
                                                             % 540     End of the header
fclose(fid);

