
function M_prop = prop_V7_vectorized_fullspec(n, d, lambda, theta2_all)
% Computes propagation matrices for all layers and wavelengths
%
% Inputs:
%   n          : 1×(N+1) refractive indices
%   d          : 1×N thicknesses
%   lambda     : 1×L wavelength vector
%   theta2_all : transmission angle(s) (scalar or 1×N)
%
% Output:
%   M_prop : 2×2×N×L propagation matrices

    nL   = n(1:end-1);         % 1×N
    d    = d(:).';             % 1×N
    Nlay = numel(nL);
    L    = numel(lambda);

    % Angles per layer
    if isscalar(theta2_all)
        th = repmat(theta2_all, 1, Nlay);  % 1×N
    else
        th = theta2_all(:).';              % 1×N
    end

    % Phase: φ_i(λ) = 2π * (n_i d_i cosθ_i) / λ
    phi_layer = (nL(:).*d(:).*cos(th(:))).';   % 1×N
    inv_lambda = 1./lambda(:).';               % 1×L
    phi = (2*pi) * (phi_layer.' * inv_lambda); % N×L

    e_neg = exp(-1i*phi);                      % N×L
    e_pos = conj(e_neg);                       % N×L

    % Pack into 2×2×N×L
    M_prop = zeros(2,2,Nlay,L,'like',phi);
    for i = 1:Nlay
        M_prop(1,1,i,:) = e_neg(i,:);
        M_prop(2,2,i,:) = e_pos(i,:);
        % off-diagonal elements remain zero
    end
end
