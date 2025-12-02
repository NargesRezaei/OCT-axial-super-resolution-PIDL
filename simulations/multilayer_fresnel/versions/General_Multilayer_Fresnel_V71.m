function [Er, Et, R, T, theta_out] = General_Multilayer_Fresnel_V71(lambda, n, d, theta_in,p, AmpSpectrum)
T(numel(lambda)) = 0;
R(numel(lambda)) = 0;
Er(numel(lambda)) = 0;
Et(numel(lambda)) = 0;
for m = 1: numel(lambda)
    [M2, theta2] = transref_V7_cell(n, theta_in,p);
    M1 = prop_V7(n, d, lambda(m), theta2);    
    M3 = cellfun(@(x,y) x*y, M2, M1, 'UniformOutput', false);
    M = eye(2);
    for h = numel(M3) : -1 : 1
        M = M * M3{h};
    end
    t = (M(1,1) * M(2,2) - M(1,2) * M(2,1)) / M(2,2);
    r = -M(2,1) / M(2,2);
    T(m) = n(end)*cos(theta2)/(n(1)*cos(theta_in))*abs(t).^2;
    R(m) = abs(r).^2;
    Et(m) = t .* AmpSpectrum(m);
    Er(m) = r .* AmpSpectrum(m);
end
theta_out = theta2;
end

%================================ SUBFUNCTIONS =====================================

function [M_cell, theta_final] = transref_V7_cell(n, theta1, p)
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
    
    % Ensure 1xN cell array shape (important for large N)
    M_cell = reshape(M_cell, 1, []);
    
    % Final transmission angle
    theta_final = theta2_all(end);
end




function M = prop_V7(n, d, lambda, theta)
% This is a function to calculate propagation matrix through a medium.
% Authors: NR
% Date: 14/03/2020
% Version 07.00
% Remove first and last element of n and convert it to a matrix with same
% size of d
%n(:, 1) = [];
n(:, end) = [];
k = 2 * pi/lambda;
phi = n .*d .* k .* cos(theta);
% Prop matrix 
M = arrayfun(@(phi_f) [exp(-1i*phi_f), 0; 0, exp(1i*phi_f)], phi, 'UniformOutput', false);
end
