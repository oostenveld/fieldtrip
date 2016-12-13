function [surface_realigned] = ft_surfacerealign(cfg, surface_original)

% FT_SURFACEREALIGN rotates, translates, scales and warps a surface mesh. The default
% is to only rotate and translate, i.e. to do a rigid body transformation in which
% only the coordinate system is changed. The different methods are described in
% detail below.
%
% INTERACTIVE - You can display the input surface mesh together with a head surface
% extracted from an anatomcial mri, and manually (using the graphical user interface)
% adjust the rotation, translation and scaling parameters, so that the surface mesh
% correspond with that of the anatomical mri.
%
% FIDUCIAL - You can apply a rigid body realignment based on three fiducial
% locations. After realigning, the fiducials in the input surface
% (typically nose, left and right ear) are along the same axes as the
% fiducials in the template surface set.
%
% HEADSHAPE - You can apply a spatial transformation/deformation that automatically
% minimizes the distance between the input surface and the head surface. The warping
% methods use a non-linear search to minimize the distance between the input surface
% vertices and their projection on the head surface.
%
% The configuration can contain the following options
%   cfg.method         = string representing the method for coregistering or aligning the electrodes
%                        'interactive'     realign manually using a graphical user interface
%                        'fiducial'        realign using three fiducials (e.g. NAS, LPA and RPA)
%   cfg.coordsys       = string specifying the origin and the axes of the coordinate
%                        system. Supported coordinate systems are 'ctf', '4d',
%                        'bti', 'yokogawa', 'asa', 'itab', 'neuromag', 'spm',
%                        'tal' and 'paxinos'. See http://tinyurl.com/ojkuhqz
%   cfg.warp          = string describing the spatial transformation for the template and headshape methods
%                        'rigidbody'       apply a rigid-body warp (default)
%                        'globalrescale'   apply a rigid-body warp with global rescaling
%                        'traditional'     apply a rigid-body warp with individual axes rescaling
%                        'nonlin1'         apply a 1st order non-linear warp
%                        'nonlin2'         apply a 2nd order non-linear warp
%                        'nonlin3'         apply a 3rd order non-linear warp
%                        'nonlin4'         apply a 4th order non-linear warp
%                        'nonlin5'         apply a 5th order non-linear warp
%
% When cfg.method = 'fiducial' and a coordinate system that is based on external
% facial anatomical landmarks (common for EEG and MEG), the following is required to
% specify the voxel indices of the fiducials:
%   cfg.fiducial.nas    = [x y z], position of Nasion
%   cfg.fiducial.lpa    = [x y z], position of LPA
%   cfg.fiducial.rpa    = [x y z], position of RPA

% Copyright (C) 2016, Simon Homoelle
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision              = '$Id$';
ft_nargin                = nargin;
ft_nargout               = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    surface_original
ft_preamble provenance surface_original
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden', 'warp');

% set the defaults
cfg.coordsys      = ft_getopt(cfg, 'coordsys');          % this allows for automatic template fiducial placement
cfg.headshape     = ft_getopt(cfg, 'headshape');         % for triangulated head surface, without labels
cfg.fiducial      = ft_getopt(cfg, 'fiducial');


if strcmp(cfg.method, 'fiducial')
  % determine the transformation matrix to align the fiducials with the axes
  transform  = ft_headcoordinates(nas, lpa, rpa, cfg.coordsys);
  % apply the transformation to the input surface
  surface_realigned = ft_warp_apply(transform, surface_original, 'homogeneous');
  % store the 4x4 transformation matrix
  surface_realigned.m = transform;
  
elseif strcmp(cfg.method, 'interactive')
  % use an external function that is shared among ft_electroderealign, ft_volumerealign and ft_surfacerealign
  tmpcfg = [];
  tmpcfg.individual.headshape = surface_original;
  tmpcfg = ft_interactiverealign(tmpcfg);
  % apply the transformation to the input surface
  surface_realigned = ft_warp_apply(tmpcfg.m, surface_original, 'homogeneous');
  % store the 4x4 transformation matrix
  surface_realigned.m = tmpcfg.m;
  
else
  error('unknown method');
end

ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   surface_original
ft_postamble provenance surface_realigned
ft_postamble history    surface_realigned
ft_postamble savevar    surface_realigned

