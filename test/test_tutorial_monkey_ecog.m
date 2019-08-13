function test_tutorial_monkey_ecog

% WALLTIME 00:20:00
% MEM 3gb

% See http://www.fieldtriptoolbox.org/tutorial/monkey_ecog/

load Event.mat
load lay
vec=1:128;
for i=1:length(vec)
 filename=strcat('ECoG_ch', num2str(vec(i)));
 data.label{i} = num2str(vec(i));
 load(filename)
 filename2 = strcat('ECoGData_ch', num2str(vec(i)));
 data.trial(i,:) = eval(filename2);
 data.time = {EventTime};
 data.fsample = 1000;
end

data.trial = {data.trial};
data.label = lay.label(1:128);
data.trial = double(data.trial{1})
data.trial = {data.trial};
data.label{129} = 'event';
clear ECoG*

%%

trigger = EventData;
sample  = EventIndex;

% determine the number of samples before and after the trigger
pretrig  = -data.fsample*2;
posttrig =  data.fsample*2;

trl = [];
for j = 2:(length(trigger)-2)
  trg1 = trigger(j);
  trg2 = trigger(j+1);
  trg3 = trigger(j+1);
  trg4 = trigger(j-1);

  %%% provided txt file reads that
  %%% orientation one spans values from 650 to 750
  %%% baseline spans values from 300 400

  if trg1 > 400 && trg2 `< 750 && trg3 >`= 650 && trg4 <= 400
    trlbegin = sample(j) + pretrig;
    trlend   = sample(j) + posttrig;
    offset   = pretrig;
    newtrl   = [trlbegin trlend offset];
    trl      = [trl; newtrl];
  end
end

%%

cfg     = [];
cfg.trl = trl;
data_1 = ft_redefinetrial(cfg,data);

%%

cfg           = [];
cfg.trl       = data_1.trial;
cfg.channel   = 'event';
cfg.ylim      = [0 3000];
cfg.blocksize = 4;
ft_databrowser(cfg, data_1);

%%

% resample the data
cfg            = [];
cfg.resamplefs = 140;
cfg.detrend    = 'no';
datads = ft_resampledata(cfg, data);

% decompose the data
cfg                 = [];
cfg.method          = 'runica';
cfg.runica.maxsteps = 50;
comp = ft_componentanalysis(cfg, datads);

% project the original data thru the components again
cfg           = [];
cfg.unmixing  = comp.unmixing;
cfg.topolabel = comp.topolabel;
comp = ft_componentanalysis(cfg, data);

%%

% perform time-frequency analysis
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:40;
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;
cfg.toi          = -2:0.05:2;
cfg.keeptrials   ='yes';
tfr = ft_freqanalysis(cfg,data);

% baseline correction
cfg               = [];
cfg.baseline      = [-1 0];
cfg.baselinetype  = 'db';
tfrbl = ft_freqbaseline(cfg, tfr);

% plot the result
cfg           = [];
cfg.channel   = {'all'};
cfg.xlim      = [-.2 2];
cfg.ylim      = [1 40];
cfg.fontsize  = 12;
cfg.layout    = lay;
% cfg.zlim    = [-.5 .5];
figure;
ft_multiplotTFR(cfg, tfrbl);

%%

% estimate high frequency gamma
cfg            = [];
cfg.output     = 'pow';
cfg.method     = 'mtmconvol';
cfg.taper      = 'dpss';
cfg.foi        = 40:2:120;
cfg.t_ftimwin  = ones(length(cfg.foi),1).*0.5;
cfg.tapsmofrq  = 6 ;
cfg.toi        = -2:0.05:2;
cfg.pad        = 'maxperlen';
cfg.keeptrials = 'yes';
tfrhf = ft_freqanalysis(cfg, data);

% baseline correct
cfg              = [];
cfg.baseline     = [-.75 0];
cfg.baselinetype = 'db';
tfrhfbl = ft_freqbaseline(cfg, tfrhf);

% plot
figure;
cfg         = [];
cfg.xlim    = [0.18 0.87]
cfg.ylim    = [53 80];
cfg.layout  = lay;
subplot(2,2,1); ft_topoplotTFR(cfg,tfrhfbl);
cfg         = [];
cfg.channel = {'chan123'};
cfg.xlim    = [-.2 1.5];
cfg.ylim    = [40 120];
subplot(2,2,2); ft_singleplotTFR(cfg,tfrhfbl);

%%

% estimate 60 Hz power
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.taper      = 'dpss';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 2;
cfg.foi        = 60;
freq           = ft_freqanalysis(cfg, data);

% then compute connectivity
cfg         = [];
cfg.method  = 'coh';
cfg.complex = 'absimag'; % check absimag solves the abs on line 161
conn = ft_connectivityanalysis(cfg,freq);

%%

% plot coherence from max gamma chan to all other
coh.label =data.label;
coh.dimord = 'chan_time'
coh.avg = conn.cohspctrm(119,:)';
coh.time = 1;

cfg           = [];
cfg.layout    = lay;
cfg.colormap  = 'jet';
cfg.zlim      = [-.2 .2];
cfg.colorbar  = 'yes';
cfg.interactive = 'no';
cfg.marker    = 'off';
cfg.highlight = 'on';
cfg.highlightchannel = {'chan123'};
cfg.highlightsymbol = '*';
cfg.highlightcolor  = [1 0 1];
cfg.highlightsize   = 12;
cfg.highlightfontsize =12;
figure;
ft_topoplotER(cfg,coh);
title('ICOH')

