function [ Zaxis, Depth ] = OCT_Analyse( OCT_Signal, Ref, Lambda, ~ )
% OCT_Analyse
% Performs OCT A-scan reconstruction from spectral interferogram.
% Steps:
%   1) Background removal
%   2) Zero-padding (spectral broadening)
%   3) λ → k linearization (resampling)
%   4) FFT to obtain depth profile

%% 1. Remove background term from interferogram
fspectrum = fftshift(fft(OCT_Signal - Ref)) / numel(OCT_Signal);

%% 2. Zero-padding in frequency domain to increase depth resolution
padorder = 4;                                   % zero-padding factor
padBord = length(fspectrum) * padorder;         % number of zeros on each side

% Broaded spectrum after adding zeros on both sides, then IFFT back
Bspectrum = ifft(ifftshift([zeros(1,padBord) fspectrum' zeros(1,padBord)])') ...
            * numel([zeros(1,padBord) fspectrum' zeros(1,padBord)]);

% New λ grid matching the zero-padded spectrum
BdLambda = linspace(Lambda(1), Lambda(end), length(Bspectrum));

%% 3. Convert wavelength → nonlinear k, then resample to linear k
K_NL = (2*pi) ./ BdLambda;                      % nonlinear wavenumber grid
K_L  = linspace(max(K_NL), min(K_NL), length(K_NL));   % linear k grid

% Interpolate spectrum onto linear-k axis
SpecInterp = interp1(K_NL, Bspectrum, K_L, 'spline');

%% 4. FFT to get depth information (A-scan)
Depth = abs(fftshift(fft(SpecInterp))) / numel(SpecInterp);

% Keep only the positive-depth half (physical part)
Depth = Depth(ceil(end/2+1):end);

%% 5. Compute depth axis
dK = K_L(1) - K_L(2);               % k sampling interval

% Maximum depth based on k-resolution
Z_max = pi / (2 * dK);

% Depth axis from 0 → Z_max
Zaxis = linspace(0, Z_max, length(Depth));

end
