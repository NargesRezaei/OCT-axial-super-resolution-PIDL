function [M_ifc, theta2_all] = transref_V7_vectorized(n, theta_in, p)
% Computes interface transfer matrices for each boundary
%
% Inputs:
%   n        : 1×(N+1) refractive indices
%   theta_in : scalar incidence angle
%   p        : 's' or 'p' polarization
%
% Outputs:
%   M_ifc      : 2×2×N interface matrices
%   theta2_all : transmission angle(s), here we return only the last one

    Nint = numel(n) - 1;

    % Snell’s law (approximation: same θ_in for all interfaces, like original code)
    theta2_all = asin(n(1:end-1).*sin(theta_in)./n(2:end));  % 1×Nint
    c1 = cos(theta_in);
    c2 = cos(theta2_all);

    n1 = n(1:end-1);
    n2 = n(2:end);

    if p == 's'   % TE
        den = (n1.*c1 + n2.*c2);
        r12 = (n1.*c1 - n2.*c2) ./ den;
        t12 = (2*n1.*c1) ./ den;

        r21 = -r12;
        t21 = (2*n2.*c2) ./ den;

    else          % 'p'  TM
        den = (n2.*c1 + n1.*c2);
        r12 = (n2.*c1 - n1.*c2) ./ den;
        t12 = (2*n1.*c1) ./ den;

        r21 = -r12;
        t21 = (2*n2.*c2) ./ (n1.*c2 + n2.*c1);
    end

    inv_t21 = 1 ./ t21;

    % Interface matrix elements
    M11 = (t12.*t21 - r12.*r21) .* inv_t21;
    M12 =  r21 .* inv_t21;
    M21 = -r12 .* inv_t21;
    M22 =  inv_t21;

    % Pack into 2×2×N array
    M_ifc = zeros(2,2,Nint,'like',M11);
    M_ifc(1,1,:) = M11;
    M_ifc(1,2,:) = M12;
    M_ifc(2,1,:) = M21;
    M_ifc(2,2,:) = M22;

    % Return last transmission angle
    theta2_all = theta2_all(end);
end
