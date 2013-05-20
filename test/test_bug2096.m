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



