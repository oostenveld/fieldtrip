function test_bug2096

% TEST test_bug2096
% TEST ft_sourcewrite




source = [];
source.dim = [5 6 7];
[X Y Z] = ndgrid(1:source.dim(1), 1:source.dim(2), 1:source.dim(3));
source.transform = eye(4);
source.pos = warp_apply(source.transform, [X(:) Y(:) Z(:)]);
source.conn = randn(prod(source.dim));
source.conndimord = 'pos_pos';

cfg = [];
cfg.filetype = 'cifti';
cfg.parameter = 'conn';
cfg.filename = 'test_bug2096';
ft_sourcewrite(cfg, source);

source = [];
source.dim = [5 6 7];
[X Y Z] = ndgrid(1:source.dim(1), 1:source.dim(2), 1:source.dim(3));
source.transform = eye(4);
source.pos = warp_apply(source.transform, [X(:) Y(:) Z(:)]);
source.time = 1:10;
source.timeseries = randn(prod(source.dim), length(source.time));
source.timeseriesdimord = 'pos_time';

cfg = [];
cfg.filetype = 'cifti';
cfg.parameter = 'timeseries';
cfg.filename = 'test_bug2096';
ft_sourcewrite(cfg, source);


% test reading these files
% http://brainvis.wustl.edu/cifti/DenseConnectome.dconn.nii
% http://brainvis.wustl.edu/cifti/DenseTimeSeries.dtseries.nii
% http://brainvis.wustl.edu/cifti/ParcellatedTimeSeries.ptseries.nii

p = '/home/common/matlab/fieldtrip/data/test/bug2096';
p = '/Volumes/Data/roboos/AeroFS/bug2096';
p = '/Users/robert/AeroFS/bug2096';
% cii1 = ft_read_mri(fullfile(p, 'DenseConnectome.dconn.nii',          'fileformat', 'cifti');
% cii2 = ft_read_mri(fullfile(p, 'DenseTimeSeries.dtseries.nii',       'fileformat', 'cifti');
% cii3 = ft_read_mri(fullfile(p, 'ParcellatedTimeSeries.ptseries.nii', 'fileformat', 'cifti');
cii1 = ft_sourceread(fullfile(p, 'DenseConnectome.dconn.nii'));
cii2 = ft_sourceread(fullfile(p, 'DenseTimeSeries.dtseries.nii'));
cii3 = ft_sourceread(fullfile(p, 'ParcellatedTimeSeries.ptseries.nii'));
cii4 = ft_sourceread(fullfile(p, 'BOLD_REST2_LR.dtseries.nii'));
cii5 = ft_sourceread(fullfile(p, 'BOLD_REST2_LR_Atlas.dtseries.nii'));


