
%  This is a function  to calculate power reflectance and transmittance for a
%  multilayer medium.
%  Authors: NR
%  Date: 03/14/2022
%  Version 02.00


function [Er, Et, R, T, theta_out] = General_Multilayer_Fresnel_V7(lambda, n, d, theta_in,p, AmpSpectrum)

T(numel(lambda)) = 0;
R(numel(lambda)) = 0;

Er(numel(lambda)) = 0;
Et(numel(lambda)) = 0;

for m = 1: numel(lambda)

    n1 = n;
    lambda1 = lambda(m);

    theta1 = theta_in;

    [M2, theta2] = transref_V7(n1, theta1,p);

    M1 = prop_V7(n1, d, lambda1, theta2);
    %M2{end+1} = eye(2);

    M3 = cellfun(@(x,y) x*y, M2, M1, 'UniformOutput',false);


    M = eye(2);
    for h = numel(M3) : -1 : 1
        M = M * M3{h};
    end


    t = (M(1,1) * M(2,2) - M(1,2) * M(2,1))/M(2,2);
    r = -M(2,1)/M(2,2);

    T(m) = n(end)*cos(theta2)/(n(1)*cos(theta_in))*abs(t).^2;
    R(m) = abs(r).^2;

     % Et(m) = T(m) * AmpSpectrum(m).^2;
     % Er(m) = R(m) * AmpSpectrum(m).^2;

    Et(m) = t .* AmpSpectrum(m);
    Er(m) = r .* AmpSpectrum(m);


end
theta_out = theta2;
end


% ================================ SUBFUNCTIONS ============================

% This is a function to calculate propagation matrix through a medium.
% Authors: NR
% Date: 14/03/2020
% Version 04.00

function M = prop_V7(n, d, lambda, theta)
% Remove first and last element of n and convert it to a matrix with same
% size of d
%n(:, 1) = [];
n(:, end) = [];


k = 2 * pi/lambda;
phi = n .*d .* k .* cos(theta);

% Prop matrix 

M = cell(1, numel(phi));
for f = 1 : numel(phi)
    M{f} = [exp(-1i.*phi(f)) 0; 0 exp(1i.*phi(f))];
end

end


% This is a function to calculate transmission and refraction matrix through amedium.
% Authors: NR
% Date: 15/03/2020
% Version 04.00

function [M, theta2] = transref_V7(n, theta1, p)

% Transfer matrix
M = cell(1, numel(n)-1);

for g = 1 : numel(n)-1
    n1 = n(g);
    n2 = n(g+1);

    theta2 = asin(n1*sin(theta1)/n2);

    if p == 's'

        r12 = (n1 * cos(theta1) - n2 * cos(theta2))/(n1 * cos(theta1) + n2 * cos(theta2));
        t12 = 2 * n1 *cos(theta1)/(n1 * cos(theta1) + n2 * cos(theta2));

        r21 = (n2 * cos(theta2) - n1 * cos(theta1))/(n2 * cos(theta2) + n1 * cos(theta1));
        t21 = 2 * n2 *cos(theta2)/(n1 * cos(theta1) + n2 * cos(theta2));
    end

    if p == 'p'
        r12 = (n2 * cos(theta1) - n1 * cos(theta2))/(n2 * cos(theta1) + n1 * cos(theta2));
        t12 = 2 * n1 *cos(theta1)/(n2 * cos(theta1) + n1 * cos(theta2));

        r21 = (n1 * cos(theta2) - n2 * cos(theta1))/(n1 * cos(theta2) + n2 * cos(theta1));
        t21 = 2 * n2 *cos(theta2)/(n1 * cos(theta2) + n2 * cos(theta1));

    end

    M{g} = 1/t21 * [t12*t21-r12*r21 r21;-r12  1];
    theta1 = theta2;
end

end
