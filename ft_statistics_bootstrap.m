function [stat, cfg] = ft_statistics_bootstrap(cfg, dat, design, varargin)

% FT_STATISTICS_BOOTSTRAP performs a nonparametric test by calculating
% bootstrap estimates. This function should not be called
% directly, instead you should call the function that is associated with the
% type of data on which you want to perform the test.
%
% Use as
%   stat = ft_timelockstatistics(cfg, data1, data2, data3, ...)
%   stat = ft_freqstatistics    (cfg, data1, data2, data3, ...)
%   stat = ft_sourcestatistics  (cfg, data1, data2, data3, ...)
%
% Where the data is obtained from FT_TIMELOCKANALYSIS, FT_FREQANALYSIS
% or FT_SOURCEANALYSIS respectively, or from FT_TIMELOCKGRANDAVERAGE,
% FT_FREQGRANDAVERAGE or FT_SOURCEGRANDAVERAGE respectively and with
% cfg.method = 'bootstrap'
%
% The configuration options that can be specified are:
%   cfg.numbootstrap     = number of bootstrap drawings
%   cfg.technique        = string, 'percentile' or 'bootstrapT' (default = XXX)
%   cfg.correctm         = string, apply multiple-comparison correction, 'no', 'max', cluster', 'bonferroni', 'holm', 'hochberg', 'fdr' (default = 'no')
%   cfg.alpha            = number, critical value for rejecting the null-hypothesis per tail (default = 0.05)
%   cfg.tail             = number, -1, 1 or 0 (default = 0)
%   cfg.feedback         = string, 'gui', 'text', 'textbar' or 'no' (default = 'text')
%   cfg.randomseed       = string, 'yes', 'no' or a number (default = 'yes')
%
% The statistic that is computed for each sample in each bootstrap drawing
% of the data is specified as
%   cfg.statistic        = 'mean'
%                          'trimmedmean'
%                          'median'
%                          'hd'          Harrell-Davis
%                          'variance'
%                          'mad'         median absolute deviation
%                          'wv'          Winsorized variance
%
% If you want to compare two or more conditions, you should specify
%   cfg.ivar             = number or list with indices, independent variable(s)
%   cfg.uvar             = number or list with indices, unit variable(s)
%
% If you use a cluster-based statistic, you can specify the following
% options that determine how the single-sample or single-voxel
% statistics will be thresholded and combined into one statistical
% value per cluster.
%   cfg.clusterstatistic = how to combine the single samples that belong to a cluster, 'maxsum', 'maxsize', 'wcm' (default = 'maxsum')
%                          option 'wcm' refers to 'weighted cluster mass',
%                          a statistic that combines cluster size and
%                          intensity; see Hayasaka & Nichols (2004) NeuroImage
%                          for details
%   cfg.clusterthreshold = method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric')
%   cfg.clusteralpha     = for either parametric or nonparametric thresholding per tail (default = 0.05)
%   cfg.clustercritval   = for parametric thresholding (default is determined by the statfun)
%   cfg.clustertail      = -1, 1 or 0 (default = 0)
%
% To include the channel dimension for clustering, you should specify
%   cfg.neighbours       = neighbourhood structure, see FT_PREPARE_NEIGHBOURS
% If you specify an empty neighbourhood structure, clustering will only be done
% in frequency and time (if available) and not over neighbouring channels.
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   cfg.correcttail      = string, correct p-values or alpha-values when doing a two-sided test, 'alpha','prob' or 'no' (default = 'no')
%   cfg.wvar             = number or list with indices, within-cell variable(s)
%   cfg.cvar             = number or list with indices, control variable(s)
%
%
%
% You can also use a custom statistic of your choise that is sensitive
% to the expected effect in the data. You can implement the statistic
% in a "statfun" that will be called for each bootstrap sample. The
% requirements on a custom statistical function is that the function
% is called statfun_xxx, and that the function returns a structure
% with a "stat" field containing the single-sample statistical values.
% Check the private functions statfun_xxx (e.g.  with xxx=tstat) for
% the correct format of the input and output.
%
% See also FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS, FT_SOURCESTATISTICS

% Undocumented local options:
%   cfg.resampling       permutation, bootstrap
%   cfg.computecritval   yes|no, for the statfun
%   cfg.computestat      yes|no, for the statfun
%   cfg.computeprob      yes|no, for the statfun
%   cfg.voxelstatistic   deprecated
%   cfg.voxelthreshold   deprecated
%   cfg.precondition     before|after|[], for the statfun

% Copyright (C) 2013, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamedval',  {'correctm', 'bonferoni', 'bonferroni'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'correctm', 'holms', 'holm'});

% set the defaults for the main function
if ~isfield(cfg, 'alpha'),               cfg.alpha    = 0.05;            end
if ~isfield(cfg, 'tail'),                cfg.tail     = 0;               end
if ~isfield(cfg, 'correctm'),            cfg.correctm = 'no';            end % no, max, cluster, bonferroni, holm, fdr
if ~isfield(cfg, 'technique'),           cfg.technique = 'percentile';   end % percentile, bootstrapT
if ~isfield(cfg, 'feedback'),            cfg.feedback = 'text';          end
if ~isfield(cfg, 'ivar'),                cfg.ivar     = 'all';           end
if ~isfield(cfg, 'uvar'),                cfg.uvar     = [];              end
if ~isfield(cfg, 'cvar'),                cfg.cvar     = [];              end
if ~isfield(cfg, 'wvar'),                cfg.wvar     = [];              end
if ~isfield(cfg, 'correcttail'),         cfg.correcttail = 'prob';       end % for the number of tails in a two-sided test
if ~isfield(cfg, 'randomseed'),          cfg.randomseed = 'yes';         end
if ~isfield(cfg, 'precondition'),        cfg.precondition = [];          end
if ~isfield(cfg, 'nicename'),            cfg.nicename = 'yes';           end

if strcmp(cfg.correctm, 'cluster')
  % set the defaults for clustering
  if ~isfield(cfg, 'clusterstatistic'),    cfg.clusterstatistic = 'maxsum';     end  % no, max, maxsize, maxsum, wcm
  if ~isfield(cfg, 'clusterthreshold'),    cfg.clusterthreshold = 'parametric'; end  % parametric, nonparametric_individual, nonparametric_common
  if ~isfield(cfg, 'clusteralpha'),        cfg.clusteralpha = 0.05;             end
  if ~isfield(cfg, 'clustercritval'),      cfg.clustercritval = [];             end
  if ~isfield(cfg, 'clustertail'),         cfg.clustertail = cfg.tail;          end
else
  % these options only apply to clustering, to ensure appropriate configs they are forbidden when _not_ clustering
  cfg = ft_checkconfig(cfg, 'unused', {'clusterstatistic', 'clusteralpha', 'clustercritval', 'clusterthreshold', 'clustertail', 'neighbours'});
end

% get the issource out flag
issource = ft_getopt(varargin, 'issource', false);

% for backward compatibility
if size(design,2)~=size(dat,2)
  design = transpose(design);
end

if ischar(cfg.ivar) && strcmp(cfg.ivar, 'all')
  cfg.ivar = 1:size(design,1);
end

% fetch function handle to the low-level statistics function
statfun = ft_getuserfun(cfg.statistic, 'statfun');
if isempty(statfun)
  error('could not locate the appropriate statistics function');
else
  fprintf('using "%s" for the single-sample statistics\n', func2str(statfun));
end

% initialize the random number generator.
if strcmp(cfg.randomseed, 'no')
  % do nothing
elseif strcmp(cfg.randomseed, 'yes')
  rand('state',sum(100*clock));
else
  % seed with the user-given value
  rand('state',cfg.randomseed);
end;

% construct the bootstrap data-shuffling matrix
fprintf('constructing randomized bootstrap samples\n');
tmpcfg = [];
tmpcfg.resampling = 'bootstrap';
tmpcfg.numresample = cfg.numbootstrap;
tmpcfg.ivar = cfg.ivar;
tmpcfg.uvar = cfg.uvar;
tmpcfg.cvar = cfg.cvar;
tmpcfg.wvar = cfg.wvar;
if isfield(cfg, 'efficient'), tmpcfg.efficient = cfg.efficient; end
resample = resampledesign(tmpcfg, design);
clear tmpcfg
Nrand = size(resample,1);

% most of the statfuns result in this warning, which is not interesting
ws = warning('off', 'MATLAB:warn_r14_stucture_assignment');

if strcmp(cfg.correctm, 'cluster')
  % determine the critical value for cluster thresholding
  if strcmp(cfg.clusterthreshold, 'nonparametric_individual') || strcmp(cfg.clusterthreshold, 'nonparametric_common')
    fprintf('using a nonparmetric threshold for clustering\n');
    cfg.clustercritval = [];  % this will be determined later
  elseif strcmp(cfg.clusterthreshold, 'parametric') && isempty(cfg.clustercritval)
    fprintf('computing a parametric threshold for clustering\n');
    tmpcfg = [];
    tmpcfg.dimord         = cfg.dimord;
    tmpcfg.dim            = cfg.dim;
    tmpcfg.alpha          = cfg.clusteralpha;
    tmpcfg.tail           = cfg.clustertail;
    tmpcfg.ivar           = cfg.ivar;
    tmpcfg.uvar           = cfg.uvar;
    tmpcfg.cvar           = cfg.cvar;
    tmpcfg.wvar           = cfg.wvar;
    if isfield(cfg, 'contrastcoefs'), tmpcfg.contrastcoefs = cfg.contrastcoefs; end % needed for Erics F-test statfun
    tmpcfg.computecritval = 'yes';  % explicitly request the computation of the crtitical value
    tmpcfg.computestat    = 'no';   % skip the computation of the statistic
    try
      cfg.clustercritval    = getfield(statfun(tmpcfg, dat, design), 'critval');
    catch
      disp(lasterr);
      error('could not determine the parametric critical value for clustering');
    end
  elseif strcmp(cfg.clusterthreshold, 'parametric') && ~isempty(cfg.clustercritval)
    fprintf('using the specified parametric threshold for clustering\n');
    cfg.clusteralpha = [];
  end
end

% compute the statistic for the observed data
ft_progress('init', cfg.feedback, 'computing statistic');
% get an estimate of the time required per evaluation of the statfun
time_pre = cputime;

try
  % the nargout function in Matlab 6.5 and older does not work on function handles
  num = nargout(statfun);
catch
  num = 1;
end

if num==1,
  % only the statistic is returned
  [statobs] = statfun(cfg, dat, design);
elseif num==2,
  % both the statistic and the (updated) configuration are returned
  [statobs, cfg] = statfun(cfg, dat, design);
else
  error('more than two statfun output variables are not supported');
end

if isstruct(statobs)
  % remember all details for later reference, continue to work with the statistic
  statfull = statobs;
  statobs  = getfield(statfull, 'stat');
else
  % remember the statistic for later reference, continue to work with the statistic
  statfull.stat = statobs;
end
time_eval = cputime - time_pre;
fprintf('estimated time per bootstrap is %.2f seconds\n', time_eval);

% pre-allocate some memory
statrand = zeros(size(statobs,1), size(resample,1));

% compute the statistic for the bootstrapped data
for i=1:Nrand
  ft_progress(i/Nrand, 'computing statistic %d from %d\n', i, Nrand);
  tmpdesign = design;                     % the design matrix is not shuffled
  tmpdat    = dat(:,resample(i,:));        % the columns of the data are resampled by means of bootstrapping
  % keep each bootstrap in memory for postprocessing
  dum = statfun(cfg, tmpdat, tmpdesign);
  if isstruct(dum)
    statrand(:,i) = getfield(dum, 'stat');
  else
    statrand(:,i) = dum;
  end
end
ft_progress('close');

switch cfg.technique
  case 'percentile'
    % FIXME this is simply counting
    prb_pos = sum(statrand<0,2);
    prb_neg = sum(statrand>0,2);
  case 'bootstrapT'
    % FIXME this is simply counting
    prb_pos = sum(statrand<0,2);
    prb_neg = sum(statrand>0,2);
  otherwise
    error('unsupported value for cfg.technique')
end

% compute the lower and upper boundaries of the confidence interval
statrand = sort(statrand,2);
stat.cilo = statrand(:, max(floor(cfg.alpha/2*Nrand), 1));
stat.cihi = statrand(:, ceil((1-cfg.alpha/2)*Nrand));

if strcmp(cfg.correctm, 'cluster')
  % do the cluster postprocessing
  [stat, cfg] = clusterstat(cfg, statrand, statobs,'issource',issource);
else
  if ~isequal(cfg.numbootstrap, 'all')
    % in case of random resampling (i.e., montecarlo sample, and NOT full
    % permutation), the minimum p-value should not be 0, but 1/N
    prb_pos = prb_pos + 1;
    prb_neg = prb_neg + 1;
    Nrand   = Nrand   + 1;
  end
  switch cfg.tail
    case 1
      clear prb_neg  % not needed any more, free some memory
      stat.prob = prb_pos./Nrand;
    case -1
      clear prb_pos  % not needed any more, free some memory
      stat.prob = prb_neg./Nrand;
    case 0
      % for each observation select the tail that corresponds with the lowest probability
      prb_neg = prb_neg./Nrand;
      prb_pos = prb_pos./Nrand;
      stat.prob = min(prb_neg, prb_pos); % this is the probability for the most unlikely tail
  end
end

% In case of a two tailed test, the type-I error rate (alpha) refers to
% both tails of the distribution, whereas the stat.prob value computed sofar
% corresponds with one tail, i.e. the probability, under the assumption of
% no effect or no difference (the null hypothesis), of obtaining a result
% equal to or more extreme than what was actually observed. The decision
% rule whether the null-hopothesis should be rejected given the observed
% probability therefore should consider alpha divided by two, to correspond
% with the probability in one of the tails (the most extreme tail). This
% is conceptually equivalent to performing a Bonferroni correction for the
% two tails.
%
% An alternative solution to distribute the alpha level over both tails is
% achieved by multiplying the probability with a factor of two, prior to
% thresholding it wich cfg.alpha.  The advantage of this solution is that
% it results in a p-value that corresponds with a parametric probability.
% Below both options are realized
if strcmp(cfg.correcttail, 'prob') && cfg.tail==0
  stat.prob = stat.prob .* 2;
  stat.prob(stat.prob>1) = 1; % clip at p=1
  % also correct the probabilities in the pos/negcluster fields
  if isfield(stat, 'posclusters')
    for i=1:length(stat.posclusters)
      stat.posclusters(i).prob = stat.posclusters(i).prob*2;
      if stat.posclusters(i).prob>1; stat.posclusters(i).prob = 1; end
    end
  end
  if isfield(stat, 'negclusters')
    for i=1:length(stat.negclusters)
      stat.negclusters(i).prob = stat.negclusters(i).prob*2;
      if stat.negclusters(i).prob>1; stat.negclusters(i).prob = 1; end
    end
  end
elseif strcmp(cfg.correcttail, 'alpha') && cfg.tail==0
  cfg.alpha = cfg.alpha / 2;
end

% compute range of confidence interval p ± 1.96(sqrt(var(p))), with var(p) = var(x/n) = p*(1-p)/N
% stddev = sqrt(stat.prob.*(1-stat.prob)/Nrand);
% stat.cirange = 1.96*stddev;

if isfield(stat, 'posclusters')
  for i=1:length(stat.posclusters)
    stat.posclusters(i).stddev  = sqrt(stat.posclusters(i).prob.*(1-stat.posclusters(i).prob)/Nrand);
    stat.posclusters(i).cirange =  1.96*stat.posclusters(i).stddev;
    if stat.posclusters(i).prob<cfg.alpha && stat.posclusters(i).prob+stat.posclusters(i).cirange>=cfg.alpha
      warning('FieldTrip:posCluster_exceeds_alpha', sprintf('The p-value confidence interval of positive cluster #%i includes %.3f - consider increasing the number of permutations!', i, cfg.alpha));
    end
  end
end
if isfield(stat, 'negclusters')
  for i=1:length(stat.negclusters)
    stat.negclusters(i).stddev  = sqrt(stat.negclusters(i).prob.*(1-stat.negclusters(i).prob)/Nrand);
    stat.negclusters(i).cirange =  1.96*stat.negclusters(i).stddev;
    if stat.negclusters(i).prob<cfg.alpha && stat.negclusters(i).prob+stat.negclusters(i).cirange>=cfg.alpha
      warning('FieldTrip:negCluster_exceeds_alpha', sprintf('The p-value confidence interval of negative cluster #%i includes %.3f - consider increasing the number of permutations!', i, cfg.alpha));
    end
  end
end

switch lower(cfg.correctm)
  case 'max'
    % the correction is implicit in the method
    fprintf('using a maximum-statistic based method for multiple comparison correction\n');
    fprintf('the returned probabilities and the thresholded mask are corrected for multiple comparisons\n');
    stat.mask = stat.prob<=cfg.alpha;
    stat.posdistribution = posdistribution;
    stat.negdistribution = negdistribution;
  case 'cluster'
    % the correction is implicit in the method
    fprintf('using a cluster-based method for multiple comparison correction\n');
    fprintf('the returned probabilities and the thresholded mask are corrected for multiple comparisons\n');
    stat.mask = stat.prob<=cfg.alpha;
  case 'bonferroni'
    fprintf('performing Bonferroni correction for multiple comparisons\n');
    fprintf('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
    stat.mask = stat.prob<=(cfg.alpha ./ numel(stat.prob));
  case 'holm'
    % test the most significatt significance probability against alpha/N, the second largest against alpha/(N-1), etc.
    fprintf('performing Holm-Bonferroni correction for multiple comparisons\n');
    fprintf('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
    [pvals, indx] = sort(stat.prob(:));                                   % this sorts the significance probabilities from smallest to largest
    k = find(pvals > (cfg.alpha ./ ((length(pvals):-1:1)')), 1, 'first'); % compare each significance probability against its individual threshold
    mask = (1:length(pvals))'<k;
    stat.mask = zeros(size(stat.prob));
    stat.mask(indx) = mask;
  case 'hochberg'
    % test the most significatt significance probability against alpha/N, the second largest against alpha/(N-1), etc.
    fprintf('performing Hochberg''s correction for multiple comparisons (this is *not* the Benjamini-Hochberg FDR procedure!)\n');
    fprintf('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
    [pvals, indx] = sort(stat.prob(:));                                   % this sorts the significance probabilities from smallest to largest
    k = find(pvals <= (cfg.alpha ./ ((length(pvals):-1:1)')), 1, 'last'); % compare each significance probability against its individual threshold
    mask = (1:length(pvals))'<=k;
    stat.mask = zeros(size(stat.prob));
    stat.mask(indx) = mask;
  case 'fdr'
    fprintf('performing FDR correction for multiple comparisons\n');
    fprintf('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
    stat.mask = fdr(stat.prob, cfg.alpha);
  otherwise
    fprintf('not performing a correction for multiple comparisons\n');
    stat.mask = stat.prob<=cfg.alpha;
end

% return the observed statistic
if ~isfield(stat, 'stat')
  stat.stat = statobs;
end

% return optional other details that were returned by the statfun
fn = fieldnames(statfull);
for i=1:length(fn)
  if ~isfield(stat, fn{i})
    stat = setfield(stat, fn{i}, getfield(statfull, fn{i}));
  end
end

if strcmp(cfg.nicename, 'yes')
  % rename the statistic and confidence interval
  stat = renamefield(stat, 'stat', cfg.statistic);
  stat = renamefield(stat, 'cilo', [cfg.statistic 'cilo']);
  stat = renamefield(stat, 'cihi', [cfg.statistic 'cihi']);
end

warning(ws); % revert to original state

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION rename a field in a structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = renamefield(s, old, new);
if isfield(s, old)
  s.(new) = s.(old);
  s = rmfield(s, old);
end
