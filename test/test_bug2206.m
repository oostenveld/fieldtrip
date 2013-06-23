% function test_bootstrap

% TEST test_bootstrap
% TEST ft_statistics_bootstrap

nrpt = 100;
nchan = 1;
ntime = 10000;

timelock1 = [];
timelock1.trial = randn(nrpt,nchan,ntime);
timelock1.label = arrayfun(@num2str, 1:nchan, 'UniformOutput', false);
timelock1.time  = 1:ntime;
timelock1.dimord = 'rpt_chan_time';

timelock2 = [];
timelock2.trial = randn(nrpt,nchan,ntime);
timelock2.label = arrayfun(@num2str, 1:nchan, 'UniformOutput', false);
timelock2.time  = 1:ntime;
timelock2.dimord = 'rpt_chan_time';

% add an effect to condition 2
effect = linspace(0,1,ntime);
for i=1:ntime
  timelock2.trial(:,:,i) = timelock2.trial(:,:,i) + effect(i);
end
timelock1.time = effect;
timelock2.time = effect;


%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.method = 'bootstrap';
cfg.statistic = 'mean';
% cfg.statistic = 'trimmedmean';
% cfg.trimpercentage = 10; % can also be specified as 0.1
cfg.numbootstrap = 1000;
bootstrap = ft_timelockstatistics(cfg, timelock1)

assert(strcmp(bootstrap.dimord), 'chan_time');

% check that the descriptive statistic is present
assert(isfield(bootstrap, 'mean'));
assert(isfield(bootstrap, 'meancilo'));
assert(isfield(bootstrap, 'meancihi'));

% these pertain to testing the difference versus zero
assert(isfield(bootstrap, 'prob'));
assert(isfield(bootstrap, 'mask'));

assert(mean(bootstrap.mask(:))>0.045 && mean(bootstrap.mask(:))<0.055); % there should be 5% false alarms 

%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.method = 'bootstrap';
cfg.technique = 'percentile'; % or bootstrapT
cfg.numbootstrap = 1000;
cfg.statistic = 'diff';
cfg.ivar = [1*ones(1,100) 2*ones(1,100)];
bootstrap = ft_timelockstatistics(cfg, timelock1, timelock2);

assert(strcmp(bootstrap.dimord), 'chan_time');

% check that the descriptive statistic is present
assert(isfield(bootstrap, 'diff'));
assert(isfield(bootstrap, 'cilo'));
assert(isfield(bootstrap, 'cihi'));

% these pertain to testing the difference between conditions
% this depends on cfg.technique
assert(isfield(bootstrap, 'prob'));
assert(isfield(bootstrap, 'mask'));


