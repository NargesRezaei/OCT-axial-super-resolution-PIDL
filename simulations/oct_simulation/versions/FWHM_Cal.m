function width = FWHM_Cal(X, K)
% Calculate FWHM of spectrum
 p = abs(abs(K) - exp(-2)*max(abs(K)));
 [m, n] = min(p);
 width = abs(X(n));
end
