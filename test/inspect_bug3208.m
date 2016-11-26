function inspect_bug3208

% TEST ft_electroderealign fT_interactiverealign ft_plot_mesh

%%

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3208'));
load headshape.mat

elec1 = ft_read_sens('standard_1020.elc');
elec2 = ft_read_sens('easycap-M10.txt');

%%

cfg = [];
cfg.elec = elec1;
cfg.method = 'interactive';
cfg.headshape = headshape;
elec = ft_electroderealign(cfg);

%%

cfg = [];
cfg.individual.elec = elec1;
cfg.template.elec = elec2;
% cfg.template.headshape = headshape;
ft_interactiverealign(cfg)

%%

cfg = [];
cfg.individual.headshape.pos = elec1.elecpos;
cfg.template.headshape = headshape;
ft_interactiverealign(cfg)


