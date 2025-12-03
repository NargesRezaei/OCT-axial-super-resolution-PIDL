%%
% OCT A-scan Simulation Program
% Author: NR
% Version: nn
% Date: 01/11/03

clc; 
clearvars; 
close all;

%% ===================== Constants =============================

theta0   = 0;
p        = 'p';
num_smpl = 150000;

%% ===================== Source Parameters =====================

lambda_c   = 840e-9;   % Center wavelength
Bandwidth  = 791e-9;   % Spectrum bandwidth
d_lambda   = 0.1e-9;   % Delta lambda
d_sbw      = 0.3e-9;
d_convert  = d_sbw / d_lambda;

lambda_i   = lambda_c - Bandwidth/2;
lambda_f   = lambda_c + Bandwidth/2;
Llambda    = Bandwidth / d_lambda + 1;
lambda     = linspace(lambda_i, lambda_f, Llambda - mod(Llambda, 2));

Llambda_sbw = Bandwidth / d_sbw + 1;
lambda_sbw  = linspace(lambda_i, lambda_f, Llambda_sbw - mod(Llambda_sbw, 2));

AmpSpectrum = gausswin(length(lambda), 39);    % Amplitude
% width_f     = 2*(lambda_c - FWHM_Cal(lambda, AmpSpectrum));
% dz          = 2*log(2)/pi * lambda_c^2 / width_f; % Axial resolution
% l_c         = dz;

%% ===================== Medium Parameters ======================

N0           = 1;
min_Ni       = 1.3;
D_N          = 0.15;
max_Ni       = min_Ni + D_N;
min_step_N   = 0.015;
N_substrate  = max_Ni;

z0           = 0;
min_OPL      = 50e-6;
D_OPL        = 500e-6;
min_step_OPL = 10e-6;

%% ===================== SBW Definition =========================

sbw_Window    = linspace(-d_sbw, (length(lambda)+2*d_convert)/10, length(lambda)+2*d_convert);
sbw_convertor = zeros(fix(numel(lambda)/d_convert), length(sbw_Window));

for kk = 1:fix(numel(lambda)/d_convert)
    sbw_convertor(kk,:) = exp(-(sbw_Window - sbw_Window(d_convert*kk + 1)).^2 ./ (2 * 0.09^2));
end

%% ===================== Spectrum Simulator =====================

Gz = fftshift(fft(AmpSpectrum)) / numel(AmpSpectrum);

%% ===================== Main Loop ==============================

% Parameters
kk = 9; % Number of layers
path = ['D:\Narges\DL-MultilayerOCT\Data\10Mm\', num2str(kk), '_3_no1stlayerlimit\'];
No_layer = kk;
N_interface = No_layer + 1;

%for jj = 1:num_smpl
% ================= Random or Custom Sample ====================
use_custom_Ni = true; % true: user-defined Ni, false: random

if use_custom_Ni
    % Example of user-defined refractive indices
    Ni_user = [1.35, 1.38, 1.2, 1.5, 1.36, 1.39, 1.34, 1.41, 1.37]; % Adjust length = No_layer
    Ni = [N0, Ni_user, N_substrate];
    OPL = linspace(20e-6, 400e-6, N_interface); % Example spacing
    OPD = [OPL(1), diff(OPL)];
    Di  = OPD ./ Ni(1:end-1);
    r12 = refl(Ni);
else
    % Random sample
    Ni(1) = min_Ni;
    OPL(1) = rand(1)*400e-6 + 15e-6;

    for g = 2:N_interface
        a = randab(min_step_OPL, (D_OPL - OPL(g-1) - (N_interface - g) * min_step_OPL));
        OPL(g) = OPL(g-1) + a;

        n2 = randab(min_Ni, max_Ni - 2*min_step_N);
        a  = (Ni(g-1) > mean([min_Ni, max_Ni]));
        b  = (abs(n2 - Ni(g-1)) < min_step_N);
        Ni(g) = n2 + b * (min_step_N - 2*a*min_step_N);
    end
    

    Ni = [N0, Ni(2:end), N_substrate];
    OPD = [OPL(1), diff(OPL)];
    Di  = OPD ./ Ni(1:end-1);

    % Shuffle layers
    s = randperm(No_layer) + 1;
    Ni(2:end-1) = Ni(s);
    Di(2:end)   = Di(s);
    OPD = Ni(1:end-1) .* Di;
    r12 = refl(Ni);

    OPL = OPD;
    for g = N_interface:-1:2
        OPL(g) = OPL(g) + sum(OPL(1:g-1));
    end
end

% ===================== OCT Simulation =====================
%{
you can use either General_Multilayer_Fresnel_V71 or General_Multilayer_Fresnel_V7. both of them have same results.
The main difference is that V7 computes each interface’s angle sequentially (layer by layer), 
while V71 uses a vectorized, one-shot calculation for all interfaces, 
which is faster but only strictly equivalent at normal incidence.
%}
[Er, Et, R, T, theta_out]        = General_Multilayer_Fresnel_V71(lambda, [N0, Ni], [z0 Di], theta0, p, AmpSpectrum/sqrt(2));
[EtR, ErR, R2, T2, theta_out2]   = General_Multilayer_Fresnel_V71(lambda, [N0, N0], z0, theta0, p, AmpSpectrum/sqrt(2));

ErR   = 0.5 .* ErR;
E_sum = Er + ErR;
I_OCT2 = smooth(0.5 * (E_sum .* conj(E_sum)));

p1     = [zeros(1,d_convert), I_OCT2', zeros(1,d_convert)];
I_OCT3 = (sbw_convertor * p1') / d_convert;

ref2   = 0.5*(Er.*conj(Er) + ErR.*conj(ErR));
p1     = [zeros(1,d_convert), ref2, zeros(1,d_convert)];
ref3   = (sbw_convertor * p1') / d_convert;

[Xaxis, Depth] = OCT_Analyse(I_OCT3, ref3, lambda_sbw , Di);
Xaxis = round(Xaxis', 7);

% ===================== Ground Truth =====================
GT = zeros(numel(Xaxis),1);
for ii = 1:numel(OPL)
    [~, h] = min(abs(Xaxis - OPL(ii)));
    GT(h, 1) = max(Gz) * r12(ii);
end

Depth = 4*Depth(1:1500)';
Xaxis = Xaxis(1:1500);

% ======================== Plot =================================
figure;
hold on;

% Plot linear scan
plot(Xaxis*1e6, Depth, 'LineWidth', 1.8);

% Plot ground truth
plot(Xaxis*1e6, GT(1:numel(Xaxis)), '.-.', 'LineWidth', 1.8);

% Labels
xlabel('Z (\mum)', 'FontSize', 14);
ylabel('Intensity (a.u.)', 'FontSize', 14);

% Legend
legend({'A-scan', 'Ground truth'}, ...
    'FontSize', 12);

grid on;
box on;
set(gca, 'FontSize', 12, 'LineWidth', 1);

% ===================== Save Data =====================
%save([path, 'ML', num2str(jj), '.mat'],'OPL', 'Depth', 'Ni');
Ni = 0; % Reset Ni
% end

% Save Xaxis
% T = table(Xaxis, 'VariableNames',{'Xaxis'});
% writetable(T,[path, 'Xaxis.txt']);
