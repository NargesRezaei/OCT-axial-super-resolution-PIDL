function [Er, Et, R, T, theta_out] = General_Multilayer_Fresnel_V11(lambda, n, d, theta_in, p, AmpSpectrum)
% GENERAL_MULTILAYER_V11 (batched over all wavelengths)
% Computes Er, Et, R, T for a multilayer stack, simultaneously for all λ
% using 3D/4D arrays and pagemtimes (or fallback loop if unavailable).
%
% Inputs:
%   lambda      : 1×L wavelength vector
%   n           : 1×(Nlay+1) refractive indices (layers + incident + exit medium)
%   d           : 1×Nlay layer thicknesses
%   theta_in    : scalar incidence angle [rad]
%   p           : polarization 's' or 'p'
%   AmpSpectrum : 1×L input amplitude spectrum
%
% Outputs (all 1×L):
%   Er, Et      : reflected / transmitted fields
%   R, T        : reflectance / transmittance
%   theta_out   : transmitted angle in the final medium (scalar)

    lambda = lambda(:).';        % force row vector
    AmpSpectrum = AmpSpectrum(:).';

    % 1) Interface matrices (2×2×Nlay) and angles
    [M_ifc, theta2] = transref_V7_vectorized(n, theta_in, p);   % 2x2xNlay

    % 2) Propagation matrices for all wavelengths (2×2×Nlay×L)
    M_prop = prop_V7_vectorized_fullspec(n, d, lambda, theta2); % 2x2xNlayxL

    % 3) Multiply interface×propagation for each layer, cascade from last to first
    L = numel(lambda);
    Nlay = size(M_ifc,3);
    M_tot = repmat(eye(2), 1, 1, L);    % 2x2xL identity matrices

    for i = Nlay:-1:1
        % M3_i = M_ifc(:,:,i) * M_prop(:,:,i,:)
        if exist('pagemtimes','builtin')
            M3_i = pagemtimes(M_ifc(:,:,i), M_prop(:,:,i,:,:));  % 2x2x1xL
            M3_i = squeeze(M3_i);                                % 2x2xL
            M_tot = pagemtimes(M_tot, M3_i);                     % 2x2xL
        else
            % Fallback for old MATLAB versions (loop over wavelengths)
            M3_i = zeros(2,2,L,'like',M_tot);
            for k = 1:L
                Mi = reshape(M_prop(:,:,i,:,k),2,2); % propagation matrix at λk
                M3_i(:,:,k) = M_ifc(:,:,i) * Mi;
                M_tot(:,:,k) = M_tot(:,:,k) * M3_i(:,:,k);
            end
        end
    end

    % 4) Extract r(λ), t(λ) from total transfer matrices
    M11 = squeeze(M_tot(1,1,:)).';
    M12 = squeeze(M_tot(1,2,:)).';
    M21 = squeeze(M_tot(2,1,:)).';
    M22 = squeeze(M_tot(2,2,:)).';

    t = (M11.*M22 - M12.*M21) ./ M22;   % 1×L
    r = -M21 ./ M22;                    % 1×L

    % 5) Compute power terms and fields
    T = (n(end)*cos(theta2)) / (n(1)*cos(theta_in)) .* abs(t).^2;  % 1×L
    R = abs(r).^2;                                                 % 1×L
    Et = t .* AmpSpectrum;                                         % 1×L
    Er = r .* AmpSpectrum;                                         % 1×L

    theta_out = theta2;   % transmitted angle in final medium
end
