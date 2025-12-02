function ErR = Reference_Mirror(lambda, N0, z0, theta0, pol, AmpSpectrum, mode, Rcoeff)
% Mirror reference at distance z0 with adjustable reflectivity
%
% Inputs:
%   lambda      - wavelength vector
%   N0          - refractive index of medium before mirror
%   z0          - distance to mirror
%   theta0      - incidence angle (radians)
%   pol         - polarization ('s' or 'p')
%   AmpSpectrum - input spectrum (field amplitude)
%   mode        - 'ideal' (default) or 'pec'
%   Rcoeff      - reflectance (0..1), e.g. 0.5 for 50% (default=1)
%
% Output:
%   ErR         - reflected reference field

if nargin < 7 || isempty(mode)
    mode = 'ideal';
end
if nargin < 8 || isempty(Rcoeff)
    Rcoeff = 1;   % default = perfect mirror
end

k = 2*pi./lambda;                                 % Wavenumber
phase = exp(-1i * 2 * k .* N0 .* z0 .* cos(theta0));  % Round-trip phase
Ein = AmpSpectrum ./ sqrt(2);                     % Split input field

% pick reflection coefficient sign depending on mode and polarization
switch lower(mode)
    case 'ideal'
        r_sign = -1;   % always -1 → π phase shift
    case 'pec'
        if pol=='s'
            r_sign = -1;   % PEC: s → -1
        else
            r_sign = +1;   % PEC: p → +1
        end
    otherwise
        error('Unknown mirror mode. Use ''ideal'' or ''pec''.');
end

% magnitude from reflectance (sqrt) and sign from mode
r = r_sign * sqrt(Rcoeff);

% reflected field
ErR = Ein .* r .* phase;
end
