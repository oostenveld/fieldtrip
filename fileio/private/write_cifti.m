function write_cifti(filename, cii)

% WRITE_CIFTI


hdr.sizeof_hdr =  540;
hdr.endian =  'b';
hdr.magic =  [110 43 50 0 13 10 26 10];
hdr.datatype =  16;
hdr.bitpix =  32;
hdr.dim =  [6 1 1 1 1 91282 1200 1];
hdr.intent_p1 =  0;
hdr.intent_p2 =  0;
hdr.intent_p3 =  0;
hdr.pixdim =  [0 1 1 1 1 1 1 1];
hdr.vox_offset =  nan;
hdr.scl_slope =  1;
hdr.scl_inter =  0;
hdr.cal_max =  0;
hdr.cal_min =  0;
hdr.slice_duration =  0;
hdr.toffset =  0;
hdr.slice_start =  0;
hdr.slice_end =  0;
hdr.descrip =  '                                                                                '; % 80
hdr.aux_file =  '                        '; % 24
hdr.qform_code =  0;
hdr.sform_code =  0;
hdr.quatern_b =  0;
hdr.quatern_c =  0;
hdr.quatern_d =  0;
hdr.qoffset_x =  0;
hdr.qoffset_y =  0;
hdr.qOffset_z =  0;
hdr.srow_x =  [0 0 0 0];
hdr.srow_y =  [0 0 0 0];
hdr.srow_z =  [0 0 0 0];
hdr.slice_code =  0;
hdr.xyzt_units =  12;
hdr.intent_code =  nan;
hdr.intent_name =  '                '; % 16
hdr.dim_info =  0;
hdr.unused_str =  '               ';  % 15

