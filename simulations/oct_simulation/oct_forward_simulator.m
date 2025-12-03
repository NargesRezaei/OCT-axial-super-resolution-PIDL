% OCT A-scan Simulation Program
clc
clearvars
close all


%% ===================== Constants ===========================
theta0       = 0;
pol          = 'p';
num_smpl     = 150000;
n_pixel_size = 1e-6;
N0           = 1.33;
z0_R         = 0;
z0_S         = 0;
Ref_scale_f  = .5;
num_smpl     = 150000;

%save_path = 'D:\Narges\DL-MultilayerOCT\Data\Retina_Simulation_V03_040601\Normal\';
%% ===================== Source Parameters =====================

% Source constants
lambda_c  = 840e-9;  %center lambda
Bandwidth = 100e-9;  %spectrum bandwidth
d_lambda  = .1e-9;   %delta lambda
FWHM      = 60e-9;

% Spectrumeter constants
d_sbw     = .3e-9;
d_convert = d_sbw/d_lambda ;

% Source Parameters
lambda_i    = lambda_c - Bandwidth/2;
lambda_f    = lambda_c + Bandwidth/2;
Llambda     = floor(Bandwidth /d_lambda)+1;
lambda      = linspace( lambda_i, lambda_f, Llambda - mod(Llambda, 2) );
Llambda_sbw = floor(Bandwidth /d_sbw)+1;
lambda_sbw  = linspace( lambda_i, lambda_f, Llambda_sbw - mod(Llambda_sbw, 2));
sigma       = FWHM / (2*sqrt(2*log(2)));
AmpSpectrum = exp( -((lambda - lambda_c).^2) / (2 * sigma^2) );
AmpSpectrum = AmpSpectrum/max(AmpSpectrum);

width_f     = 2*(lambda_c - FWHM_Cal(lambda, AmpSpectrum));
dz          = 2 * log(2)/pi* lambda_c^2/width_f; % axial resolution of OCT
l_c         = dz;


%% ======================== Medium Parameter ====================
% N
N0           = 1;
min_Ni       = 1.3;
D_N          = 0.15;
max_Ni       = min_Ni + D_N;
min_step_N   = 0.015;
N_substrate  = max_Ni;

% OPL
z0           = 0;
min_OPL      = 50e-6;
D_OPL        = 500e-6;
min_step_OPL = 10e-6;

%% ======================== Define SBW ==========================

sbw_Window = linspace(lambda_i, lambda_f, Llambda_sbw);
sbw_convertor = zeros(Llambda_sbw, length(lambda));
for kk = 1 : Llambda_sbw
    sbw_convertor(kk,:) = exp(-((lambda - lambda_sbw(kk)).^2) / (2 * (d_sbw/2)^2));
    sbw_convertor(kk,:) = sbw_convertor(kk,:) / sum(sbw_convertor(kk,:));
end

Gz = fftshift(fft(AmpSpectrum))/(numel(AmpSpectrum));

%% ======================== Reference Arm ==========================
ErR = Reference_Mirror(lambda, N0, z0_R, theta0, pol, AmpSpectrum, 'ideal', Ref_scale_f);

%% ======================== Sample Arm =============================


for kk = 9
    path        = ['D:\Narges\DL-MultilayerOCT\Data\10Mm\', num2str(kk), '_3_no1stlayerlimit\'];
    No_layer    = kk;
    N_interface = No_layer + 1;

    %%for jj = 1 : num_smpl

    OPD = [1 1];
    r12 = [1 1];

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

        Ni(1) = min_Ni;
        OPL(1) = rand(1)*400e-6 + 15e-6;
        for g = 2 : N_interface
            a      = randab(min_step_OPL, (D_OPL - OPL(g-1) - (N_interface - g) * min_step_OPL));
            OPL(g) = OPL(g-1) + a;

            n2 = randab(min_Ni, max_Ni - 2*min_step_N);
            a  = (Ni(g-1) > mean([min_Ni, max_Ni]));
            b  = (abs(n2 - Ni(g-1)) < min_step_N);
            Ni(g) = n2 + b * (min_step_N - 2*a*min_step_N);
        end

        Ni          = [N0 Ni(2:end), N_substrate];
        OPD         = [OPL(1) OPL(2:end) - OPL(1:end-1)];
        Di          = OPD ./ Ni(1:end-1);
        s           = randperm(No_layer)+1;
        Ni(2:end-1) = Ni(s);
        Di(2:end)   = Di(s);
        OPD         = Ni(1:end-1).*Di;
        r12         = refl(Ni);

        OPL = OPD;
        for g = N_interface :-1: 2
            OPL(g) = OPL(g) + sum(OPL(1:g-1));
        end
    end

    %% ======================== Spectral Interferogram =============================
    [Er, ~, ~, ~, ~] = General_Multilayer_V11(lambda, [N0, Ni],  [z0_S  Di], theta0, pol, AmpSpectrum/sqrt(2));
    E_sum = Er +  ErR;
    I_OCT2 = 1/2 * (E_sum .* conj(E_sum));

    %% ======================== Spectrometer =======================================
    I_OCT3 = sbw_convertor * I_OCT2';

    ref2 = 1/2 * (Er .* conj(Er) + ErR .* conj(ErR));
    ref3 = sbw_convertor * ref2';

    %% ======================== Spectral Interferogram Analyser ====================

    [ Xaxis, Depth ] = OCT_Analyse(I_OCT3, ref3, lambda_sbw , Di);
    Xaxis = round(Xaxis', 7);



    GT = zeros(numel(Xaxis),1);
    for ii = 1 : numel(OPL)
        h = find(abs(Xaxis - OPL(ii)) == min(abs(Xaxis - OPL(ii))));
        GT(h, 1) = max(Gz)*r12(ii)/4;
    end


    Depth = 4*Depth(1:200)';
    Xaxis = Xaxis(1:200);
    %% ======================== Plot =============================================
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


    %% ======================== Save =============================================
    % save([path, 'ML', num2str(jj), '.mat'], 'OPL', 'Depth', 'Ni');

    Ni = 0; % Reset Ni

end
% end

% Save Xaxis
% T = table(Xaxis, 'VariableNames',{'Xaxis'});
% writetable(T,[path, 'Xaxis.txt']);
