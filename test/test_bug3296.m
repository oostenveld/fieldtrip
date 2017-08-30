function test_bug2978

% WALLTIME 00:20:00
% MEM 2GB

% TEST ft_multiplotER


dataset = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare the data

cfg = [];
cfg.channel = 'MEG';
cfg.demean = 'yes';
cfg.dataset = dataset;
data = ft_preprocessing(cfg); % 266 trials

cfg = [];
cfg.trials = 1:10;
data = ft_selectdata(cfg, data); % 10 trials to speed it all up
data.trial{10}(:,450:end) = data.trial{10}(:,450:end) + 1e-11; % insert a jump in trial 10

cfg = [];
timelock1 = ft_timelockanalysis(cfg, data);
cfg.keeptrials = 'yes';
timelock2 = ft_timelockanalysis(cfg, data);

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
freq1 = ft_freqanalysis(cfg, data);
cfg.keeptrials = 'yes';
freq2 = ft_freqanalysis(cfg, data);

cfg = [];
cfg.foilim = [1 20];
cfg.method = 'wavelet';
cfg.toi = -1:0.100:2;
timefreq1 = ft_freqanalysis(cfg, data);
cfg.keeptrials = 'yes';
timefreq2 = ft_freqanalysis(cfg, data);

% I don't know how to get single trial (pseudo) estimates of coherence
if false
  cfg = [];
  cfg.method = 'mtmfft';
  cfg.taper = 'hanning';
  cfg.output = 'fourier'; % this implies keeptrials
  powandcsd = ft_freqanalysis(cfg, data);
  
  cfg = [];
  cfg.method = 'coh';
  coh1 = ft_connectivityanalysis(cfg, powandcsd);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the figures, timelock

cfg = [];
cfg.layout = 'CTF151.lay';

cfg.trials = 1;
figure; ft_multiplotER(cfg, timelock1);
cfg.trials = 'all';
figure; ft_multiplotER(cfg, timelock1); % same as previous
cfg.trials = 10;
figure; ft_multiplotER(cfg, timelock1); % error
cfg.trials = [];
figure; ft_multiplotER(cfg, timelock1); % error

cfg.trials = 1;
figure; ft_multiplotER(cfg, timelock2); % no jump
cfg.trials = 'all';
figure; ft_multiplotER(cfg, timelock2); % small jump (average)
cfg.trials = 10;
figure; ft_multiplotER(cfg, timelock2); % small jump
cfg.trials = [];
figure; ft_multiplotER(cfg, timelock2); % error

cfg.trials = 1;
figure; ft_multiplotER(cfg, data); % no jump
cfg.trials = 'all';
figure; ft_multiplotER(cfg, data); % small jump (average)
cfg.trials = 10;
figure; ft_multiplotER(cfg, data); % big jump
cfg.trials = [];
figure; ft_multiplotER(cfg, data); % error


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the figures, power

cfg = [];
cfg.layout = 'CTF151.lay';

cfg.trials = 1;
figure; ft_multiplotER(cfg, freq1);
cfg.trials = 'all';
figure; ft_multiplotER(cfg, freq1);
cfg.trials = 10;
figure; ft_multiplotER(cfg, freq1);
cfg.trials = [];
figure; ft_multiplotER(cfg, freq1);

cfg.trials = 1;
figure; ft_multiplotER(cfg, freq2);
cfg.trials = 'all';
figure; ft_multiplotER(cfg, freq2);
cfg.trials = 10;
figure; ft_multiplotER(cfg, freq2);
cfg.trials = [];
figure; ft_multiplotER(cfg, freq2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the figures, timefreq

cfg = [];
cfg.layout = 'CTF151.lay';

cfg.trials = 1;
figure; ft_multiplotTFR(cfg, timefreq1);
cfg.trials = 'all';
figure; ft_multiplotTFR(cfg, timefreq1);
cfg.trials = 10;
figure; ft_multiplotTFR(cfg, timefreq1);
cfg.trials = [];
figure; ft_multiplotTFR(cfg, timefreq1);

cfg.trials = 1;
figure; ft_multiplotTFR(cfg, timefreq2);
cfg.trials = 'all';
figure; ft_multiplotTFR(cfg, timefreq2);
cfg.trials = 10;
figure; ft_multiplotTFR(cfg, timefreq2);
cfg.trials = [];
figure; ft_multiplotTFR(cfg, timefreq2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the figures, bivariate

if false
  % FIXME construct appropriate bivariate data
  cfg = [];
  cfg.layout = 'CTF151.lay';
  cfg.parameter = 'cohspctrm';
  cfg.refchannel = 42;
  figure; ft_multiplotER(cfg, coh1);
end