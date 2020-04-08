function [glm] = ft_glmanalysis(cfg, data, model)

% FT_GLMANALYSIS fits a General Linear Model to the data
%
% Use as
%   glm = ft_glmanalysis(cfg, data)
% or as
%   glm = ft_glmanalysis(cfg, data, model)
% where the input data and the optional model are both formatted a raw data
% structure, similar to the output of FT_PREPROCESSING.
%
% Specifying only "data" as a single input argument makes sense if it includes
% channels that describe your hypothetical model, or if it includes reference
% channels that you want to regress out of your data. Specifying both input "data"
% and "model" makes sense if your model underwent different processing than your
% data, and/or if your model contains channels that are also present in your data
% itself.
%
% The configuration should contain
%
%   cfg.channel     = Nx1 cell-array with selection of channels (default = 'all'),
%                      see FT_CHANNELSELECTION for details
%   cfg.glmchannel  = the channels with the model specification (default is all channels from the "model" input)
%   cfg.normalize   = string, 'yes' or 'no', normalization to make the confounds
%                     orthogonal (default = 'yes')
%   cfg.output      = 'residual' (default), 'beta', or 'model'.
%                     If 'residual' is specified, the output is a data
%                     structure containing the residuals after regressing
%                     out the in cfg.reject listed confounds. If 'beta' or 'model'
%                     is specified, the output is a data structure containing
%                     the regression weights or the model, respectively.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_REGRESSCONFOUND, FT_DENOISE_PCA, FT_TIMELOCKSTATISTICS

% Copyright (C) 2020, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function

ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    data model
ft_preamble provenance data model
ft_preamble trackconfig

if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');
if ft_nargin>2
  model = ft_checkdata(model, 'datatype', {'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');
end

% specify the defaults
cfg.trials     = ft_getopt(cfg, 'trials',      'all', 1);
cfg.channel    = ft_getopt(cfg, 'channel',     'all');
cfg.normalize  = ft_getopt(cfg, 'normalize', 'yes');
cfg.output     = ft_getopt(cfg, 'output', 'residual');
cfg.statistics = ft_getopt(cfg, 'statistics', 'no');
cfg.ftest      = ft_getopt(cfg, 'ftest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize the data

if ft_nargin>2
  % select trials and channels/regressors from the second input argument
  tmpcfg = keepfields(cfg, {'trials'});
  tmpcfg.channel = cfg.glmchannel; % this one is named differently
  model = ft_selectdata(tmpcfg, model);
  % do not restore the provenance information
else
  % select trials and channels/regressors from the first input argument
  tmpcfg = keepfields(cfg, {'trials'});
  tmpcfg.channel = cfg.glmchannel; % this one is named differently
  model = ft_selectdata(tmpcfg, data);
  % do not restore the provenance information
end

% organize the regressors in a Nsamples*Nchans matrix
regr = cat(2, model.trial{:})';

% select trials and channels of interest, this must be done after getting the regressors out
tmpcfg = keepfields(cfg, {'trials', 'channel'});
data   = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

% organize the data in a Nsamples*Nchans matrix
dat = cat(2, data.trial{:})';

% regressor normalization for orthogonality
if strcmp(cfg.normalize, 'yes')
  fprintf('normalizing the regressors, except the constant \n');
  for c = 1:size(regr,2)
    AVG = mean(regr(:,c));
    STD = std(regr(:,c),0,1);
    if abs(STD/AVG)<10*eps
      fprintf('confound %s is a constant \n', num2str(c));
    else
      regr(:,c) = (regr(:,c) - AVG) / STD;
    end
    clear AVG STD;
  end
else
  fprintf('skipping normalization procedure \n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLM MODEL
%   Y = X * B + err, where Y is data, X is the model, and B are beta's
% which means
%   Best = X\Y ('matrix division', which is similar to B = inv(X)*Y)
% or when presented differently
%   Yest = X * Best
%   Yest = X * X\Y
%   Yclean = Y - Yest (the true 'clean' data is the recorded data 'Y' -
%   the data containing confounds 'Yest')
%   Yclean = Y - X * X\Y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% estimate and remove the confounds
fprintf('estimating the regression weights\n');

if ~any(isnan(dat(:))) % if there are no NaNs, process all at once
  beta = regr\dat;                                                        % B = X\Y
else % otherwise process per colum set as defined by the nan distribution
  [u,i,j] = unique(~isnan(dat)','rows','first');      % find unique rows
  uniquecolumns = u';                                 % unique column types
  Nuniques      = numel(i);                           % number of unique types
  beta_temp     = NaN(Nuniques, nconf, size(dat,2));  % declare empty variable
  for n = 1:Nuniques                        % for each unique type
    rowidx = find(uniquecolumns(:,n)==1);   % row indices for unique type
    colidx = find(j==n);                    % column indices for unique type
    if any(uniquecolumns(:,n))              % if vector contains a nonzero number
      beta_temp(n,:,colidx) = regr(rowidx,:)\dat(rowidx,colidx);         % B = X\Y
    end
  end
  beta = reshape(nansum(beta_temp,1),[nconf size(dat,2)]); % sum the betas
  clear beta_temp
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize the output
glm = [];

switch cfg.output
  case 'beta'
    % the output is not a standard FieldTrip structure
    glm.beta        = beta';
    glm.betadimord  = 'chan_regressor';
    glm.label       = data.label;
    glm.regressor   = model.label;
    
  case 'comp'
    % the output is according to FT_DATATYPE_COMP
    glm.topo      = beta';
    glm.topolabel = data.label;
    glm.trial     = model.trial;
    glm.label     = model.label;

  case 'residual'
    % the output is according to FT_DATATYPE_RAW
    glm.time  = data.time;
    glm.label = data.label;

    endsample = nan(1, numel(data.trial));
    for i=1:numel(data.trial)
      endsample(i) = size(data.trial{i},2);
    end
    endsample = cumsum(endsample);
    begsample = [1 endsample(1:end-1)-1];
    
    residual = (dat - regr*beta)';
    for i=1:numel(data.trial)
      glm.trial{i} = residual(:,begsample(i):endsample(i));
    end
    
  case 'model'
    % the output is according to FT_DATATYPE_RAW
    glm.time  = data.time;
    glm.label = data.label;

    endsample = nan(1, numel(data.trial));
    for i=1:numel(data.trial)
      endsample(i) = size(data.trial{i},2);
    end
    endsample = cumsum(endsample);
    begsample = [1 endsample(1:end-1)-1];
    
    model = (regr*beta)';
    for i=1:numel(data.trial)
      glm.trial{i} = model(:,begsample(i):endsample(i));
    end
    
  otherwise
    error('output ''%s'' is not supported', cfg.output);
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data model
ft_postamble provenance glm
ft_postamble history    glm
ft_postamble savevar    glm
