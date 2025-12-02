function [Er, Et, R, T, theta_out] = General_Multilayer_Fresnel_V10(lambda, n, d, theta_in, p, AmpSpectrum)
% Calculates power reflectance (R), transmittance (T), reflected electric field (Er),
% transmitted electric field (Et), and output angle (theta_out) for a
% multilayer structure., 
%
% Inputs:
%   lambda: Wavelengths [vector]
%   n: Refractive indices [vector]
%   d: Layer thicknesses [vector]
%   theta_in: Incidence angle [radians]
%   p: Polarization ('s' or 'p')
%   AmpSpectrum: Input electric field amplitude spectrum
%
% Outputs:
%   Er: Reflected field (complex)
%   Et: Transmitted field (complex)
%   R: Power reflectance
%   T: Power transmittance
%   theta_out: Output transmission angle

% Initialize output arrays
T = zeros(size(lambda));
R = zeros(size(lambda));
Er = zeros(size(lambda));
Et = zeros(size(lambda));

q = numel(lambda);
% Step 1: Calculate interface matrices (M2) and propagation angles (theta2)
[M2, theta2] = transref_V7_vectorized(n, theta_in, p, q);

% Step 2: Calculate propagation matrices (M1) for each layer
M1 = prop_V7_vectorized(n, d, lambda, theta2);

% Step 3: Multiply M2 and M1 matrices element-wise
M3 = cellfun(@(x,y) x*y, M2, M1, 'UniformOutput', false);
M3 = M3';
M_M = reshape(M3, 1, []);  % Flatten to 1xN cell array

% Step 4: Cascade matrix multiplication with downsampling

M = eye(2);
M_f{q} = M;
for k = numel(M_M):-1:1
    M = M * M_M{k};
    if mod(k, size(M3, 1)) == 1 || size(M3, 1) == 1 % Downsample at original layer boundaries
        M_f{q} = M;
        q = q - 1;
        M = eye(2);  % Reset for next block
    end
end

% Step 5: Calculate transmission (t) and reflection (r) coefficients
t = cell2mat(cellfun(@(M) (M(1,1)*M(2,2) - M(1,2)*M(2,1))/M(2,2), M_f, 'UniformOutput', false));
r = cell2mat(cellfun(@(M) -M(2,1)/M(2,2), M_f, 'UniformOutput', false));

% Step 6: Compute power terms and output fields
T = n(end)*cos(theta2)/(n(1)*cos(theta_in)) .* abs(t).^2;  % Transmittance
R = abs(r).^2;  % Reflectance
Et = t .* AmpSpectrum;  % Transmitted field
Er = r .* AmpSpectrum;  % Reflected field

theta_out = theta2;  % Output angle
end

%======================== SUBFUNCTIONS ========================%

function [M_cell, theta_final] = transref_V7_vectorized(n, theta1, p, q)
    % VECTORIZED VERSION OF TRANSREF_V7
    % Returns a 1xN cell array where each cell contains a 2x2 transfer matrix
    % Inputs:
    %   n - refractive indices vector
    %   theta1 - initial angle of incidence (radians)
    %   p - polarization ('s' or 'p')
    % Outputs:
    %   M_cell - 1x(N-1) cell array of 2x2 matrices
    %   theta2 - final transmission angle
    
    % Calculate for all layers simultaneously
    g_indices = 1:numel(n)-1;  % Indices for layer interfaces
    
    % Vectorized Snell's law calculation
    theta2_all = asin(n(g_indices).*sin(theta1)./n(g_indices+1));
    
    % Precompute cosine terms for efficiency
    cos_theta1 = cos(theta1);
    cos_theta2_all = cos(theta2_all);
    
    % Polarization-specific calculations
    if p == 's'
        % TE (s-polarization) Fresnel coefficients
        r12 = (n(g_indices).*cos_theta1 - n(g_indices+1).*cos_theta2_all) ./ ...
              (n(g_indices).*cos_theta1 + n(g_indices+1).*cos_theta2_all);
        t12 = (2 * n(g_indices).*cos_theta1) ./ ...
              (n(g_indices).*cos_theta1 + n(g_indices+1).*cos_theta2_all);
        r21 = -r12;  % Simplified relation for s-pol
        t21 = (2 * n(g_indices+1).*cos_theta2_all) ./ ...
              (n(g_indices).*cos_theta1 + n(g_indices+1).*cos_theta2_all);
    elseif p == 'p'
        % TM (p-polarization) Fresnel coefficients
        r12 = (n(g_indices+1).*cos_theta1 - n(g_indices).*cos_theta2_all) ./ ...
              (n(g_indices+1).*cos_theta1 + n(g_indices).*cos_theta2_all);
        t12 = (2 * n(g_indices).*cos_theta1) ./ ...
              (n(g_indices+1).*cos_theta1 + n(g_indices).*cos_theta2_all);
        r21 = -r12;  % Simplified relation for p-pol
        t21 = (2 * n(g_indices+1).*cos_theta2_all) ./ ...
              (n(g_indices).*cos_theta2_all + n(g_indices+1).*cos_theta1);
    end
    
    % Vectorized transfer matrix calculation
    % Each matrix M = (1/t21)*[t12*t21-r12*r21, r21; -r12, 1]
    M_cell = arrayfun(@(r12,r21,t12,t21) (1/t21)*[t12*t21-r12*r21, r21; -r12, 1], ...
                     r12, r21, t12, t21, 'UniformOutput', false);

    M_cell = repmat(M_cell, [q, 1]);
    % Ensure 1xN cell array shape (important for large N)
    %M_cell = reshape(M_cell, 1, []);
    
    % Final transmission angle
    theta_final = theta2_all(end);
end

function M = prop_V7_vectorized(n, d, lambda, theta)
% PROPAGATION MATRIX CALCULATOR
% Computes 2x2 propagation matrices for all layers/wavelengths
%
% Inputs:
%   n: Refractive indices (excluding last layer)
%   d: Layer thicknesses
%   lambda: Wavelengths [vector]
%   theta: Propagation angles [vector]
%
% Output:
%   M: Cell array of propagation matrices

% Remove last refractive index (not needed for propagation)
n(:, end) = [];

% Wavenumber
k = 2*pi./lambda;

% Phase accumulation
phi = n .* d .* cos(theta) .* k';

% Create diagonal propagation matrices
M = arrayfun(@(phi_f) [exp(-1i*phi_f), 0; 0, exp(1i*phi_f)], phi, 'UniformOutput', false);
end
