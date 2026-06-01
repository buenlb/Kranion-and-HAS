function Params = setInsightecERFA_Params(varargin)
% FUSF_ERFAParams  Returns ERFA parameters for the Insightec ExAblate 4000.
%                  All parameters have defaults but can be overridden using
%                  name-value pairs.
%
% USAGE:
%   Params = FUSF_ERFAParams()
%   Params = FUSF_ERFAParams('fMHz', 0.72, 'R', 0.14)

% Default Parameters
% Transducer and ERFA parameters for Insightec ExAblate
Params.fMHz   = 0.6700;          % Frequency (MHz)
Params.R      = 0.1500;          % Radius of curvature (m)
Params.Dv     = 0.22;            % Arc dimension, vertical (m)
Params.Dh     = 0.22;            % Arc dimension, horizontal (m)
Params.imax   = 301;             % Number of discretization points along Dv (Less than 1% change if increased to 401)
Params.kmax   = 301;             % Number of discretization points along Dh (Less than 1% change if increased to 401)
Params.relem  = 0.0055;          % Element radius (m)
Params.d      = 0.0532;          % Distance from transducer back to ERFA plane (m)
Params.Lv     = 0.25;            % ERFA plane size, vertical (m)
Params.Lh     = 0.25;            % ERFA plane size, horizontal (m)
Params.lmax   = 501;             % ERFA plane increments, vertical (must be odd)
Params.mmax   = 501;             % ERFA plane increments, horizontal (must be odd)
Params.isPA   = 1;               % Phased array flag (1=PA, 0=solid)
Params.c0     = 1500;            % Speed of sound in water (m/s)
Params.rho0   = 1000;            % Density of water (kg/m^3)
Params.sName  = 'ERFA8.mat';     % Full save name for saving ouput

if mod(length(varargin),2)
    error('Must be in parameter, value pairs')
end

for i = 1:2:length(varargin)
    if isfield(Params, varargin{i})
        Params.(varargin{i}) = varargin{i+1};
    else
        error('Unknown parameter: %s', varargin{i});
    end
end

if Params.imax~=Params.kmax
    warning('Different discretization in theta and phi.')
end

if Params.lmax~=Params.mmax
    warning('ERFA is not square!')
end