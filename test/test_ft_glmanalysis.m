function test_ft_glmanalysis

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_glmanalysis

%%

data = [];
data.fsample = 1000;
for i=1:10
  data.label{i} = sprintf('chan%02i', i);
end
for i=1:10
  data.time{i}  = (0:999)/1000;
  data.trial{i} = randn(10, 1000);
end
data.label = data.label(:);

model = [];
model.fsample = 1000;
for i=1:5
  model.label{i} = sprintf('model%02i', i);
end
for i=1:10
  model.time{i}   = (0:999)/1000;
  model.trial{i} = randn(5, 1000);
end
model.label = model.label(:);

both = ft_appenddata([], data, model);


%%

cfg = [];
cfg.channel = data.label;
cfg.glmchannel = model.label;
glm1a = ft_glmanalysis(cfg, data, model);

glm1b = ft_glmanalysis(cfg, both); 

% should have the same result, only the cfg will be slightly different
assert(isequal(rmfield(glm1a, 'cfg'), rmfield(glm1b, 'cfg')))

%%

cfg = [];
cfg.channel = data.label;
cfg.glmchannel = model.label;

cfg.output = 'residual';
residual = ft_glmanalysis(cfg, data, model);

cfg.output = 'model';
estimated = ft_glmanalysis(cfg, data, model);

cfg.output = 'beta';
beta = ft_glmanalysis(cfg, data, model);

cfg.output = 'comp';
comp = ft_glmanalysis(cfg, data, model);
