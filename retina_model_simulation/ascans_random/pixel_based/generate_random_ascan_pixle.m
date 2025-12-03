% Genrating Retina A-scan
% Athour: NR
% V 00
% 04/03/30

clc
clearvars
close all


save_path = 'D:\Narges\DL-MultilayerOCT\Data\Retina_Simulation_V02_040426\Abnormal\';
%%%%%%%%%%%%%%%%%%%%
%%%       constants
%%%%%%%%%%%%%%%%%%%%

theta0       = 0;
pol          = 'p';
num_smpl     = 100000;
pixel_size   = 2.8e-6;        % pixel size in micrometers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% medium parameters

% N & OPL
N0           = 1.33;
z0           = 0;
N_1          = [N0 1.38, 1.37, 1.376, 1.39, 1.398, 1.41, 1.425, 1.405, 1.395, 1.385 ,  N0];
Depth_ranges = [10 15; 3 60; 3 60; 3 60; 3 60; 3 60; 3 10; 3 60; 3 10; 3 60; 3 60]*1e-6;
dn           = [0 0.1 0.02 0.005 0.02 0.003 0.01 0.003 0.2 0.008 0.2 0];         % maximum deviation for refractive index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spectrume parameters

% Source constants
lambda_c  = 840e-9;  %center lambda
Bandwidth = 100e-9;  %spectrum bandwidth
d_lambda  = .1e-9;   %delta lambda
FWHM = 60e-9;

% Spectrumeter constants
d_sbw     = .3e-9;
d_convert = d_sbw/d_lambda ;

% Source Parameters
lambda_i    = lambda_c - Bandwidth/2;
lambda_f    = lambda_c + Bandwidth/2;
Llambda     = Bandwidth /d_lambda+1;
lambda      = linspace( lambda_i, lambda_f, Llambda - mod(Llambda, 2) );
Llambda_sbw = Bandwidth /d_sbw+1;
lambda_sbw  = linspace( lambda_i, lambda_f, Llambda_sbw - mod(Llambda_sbw, 2));
sigma       = FWHM / (2*sqrt(2*log(2)));
AmpSpectrum = exp( -((lambda - lambda_c).^2) / (2 * sigma^2) );
AmpSpectrum = AmpSpectrum/max(AmpSpectrum);

width_f     = 2*(lambda_c - test02(lambda, AmpSpectrum));
dz          = 2 * log(2)/pi* lambda_c^2/width_f; % axial resolution of OCT
l_c         = dz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Define SBW

sbw_Window    = linspace(-d_sbw, (length(lambda)+2*d_convert)/10, length(lambda)+2*d_convert);
sbw_convertor = zeros(fix((numel(lambda))/d_convert), length(sbw_Window));
for kk = 1 : fix((numel(lambda))/d_convert)
    sbw_convertor(kk,:) = exp(-(sbw_Window-sbw_Window(d_convert*kk+1)).^2./(2 * 0.09^2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reference arm
[EtR, ErR, R2, T2, theta_out2] = General_Multilayer_V10(lambda, [N0, N0], z0, theta0, pol, AmpSpectrum/sqrt(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sample arm
pxl_ranges = round(Depth_ranges / pixel_size);

for jj = 1 : num_smpl
    % Pixel range (converted from micrometers to pixels)
    num_layers = randsample([9, 10, 11], 1, true, [0.1, 0.1, 0.8]);
    pxl_s = pxl_ranges(end-num_layers+1:end,:);        % number of layers

    % number of layers
    N_layer     = [N_1(1) N_1(end-num_layers+1:end)];
    dn_layer    = [dn(1) dn(end-num_layers+1:end)];

    % Initialize cells for thicknesses and refractive indices
    thick_cells  = cell(1, num_layers);
    refidx_cells = cell(1, num_layers);
    OPL          = zeros(1,size(num_layers, 1));


    len1 = randi(pxl_s(1, :));
    thick_cells{1}  = pixel_size * ones(1, len1);
    refidx_cells{1} = N_layer(1) * ones(1, len1);
    OPL(1)          = sum(thick_cells{1} .* refidx_cells{1});
    N_trans = 3;
    for i = 2:num_layers
        
        len = randi(pxl_s(i, :));

        base_value = N_layer(i);
        ref_index_raw = (2 * rand(1, len) - 1) * dn_layer(i)/2 + base_value;
        
        % Check if len < N_trans
        num_trans = min(N_trans, len);
        
        if i > 3 && i<11 
            prev_base = N_layer(i-1);
            
            % Transition zone
            for k = 1:num_trans
                alpha = k / num_trans;
                ref_index_raw(k) = prev_base*(1-alpha) + base_value*alpha + ...
                    (2*rand - 1)*dn_layer(i)/4;
            end
        end

        % Smoothing
        ref_index = ref_index_raw;

        % Thickness
        thick_array = pixel_size * ones(1,len);

        % Ensure same length
        % if length(ref_index) ~= length(thick_array)
        %     error('Length mismatch in layer %d: ref_index=%d, thick_array=%d', ...
        %         i,length(ref_index),length(thick_array));
        % end

        refidx_cells{i} = ref_index;
        thick_cells{i}  = thick_array;

        OPL(i) = sum(thick_array .* ref_index);
    end


    % Concatenate all thickness and refractive index values
    LD = cell2mat(thick_cells);         % total layer thickness array
    LN = cell2mat(refidx_cells);       % total refractive index array

    % Append final refractive index value (bottom layer)
    LN = [LN, N_1(end)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spectral Interferogram
    [Er, Et, R, T, theta_out] = General_Multilayer_V10(lambda, [N0, LN],  [z0  LD], theta0, pol, AmpSpectrum/sqrt(2));
    ErR_2                     = .5 .* ErR;
    E_sum                     = Er +  ErR_2;
    I_OCT2                    = smooth(1/2*(E_sum .* conj(E_sum)));
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spectrometer 
    p1      = [zeros(1,d_convert) I_OCT2' zeros(1,d_convert)];
    I_OCT3  = (sbw_convertor*p1')/d_convert;

    ref2  =  1/2*(Er .* conj(Er) +  ErR_2 .* conj(ErR_2));
    p1    = [zeros(1,d_convert) ref2 zeros(1,d_convert)];
    ref3  =  (sbw_convertor*p1')/d_convert;

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spectral Interferogram Analyse
    % [ Xaxis, Depth ] = OCT_Analyse(I_OCT2, ref2', lambda , LD);
    [ Xaxis, Depth ] = OCT_Analyse(I_OCT3, ref3, lambda_sbw , LD);

    Xaxis = round(Xaxis', 7);

    Depth = Depth(1:200)';
    Xaxis = Xaxis(1:200);
    OPL_2 = cumsum(OPL);
   
        figure
    cla
    plot(Xaxis*1e6, Depth)
    xlabel('Z (\mum)', 'FontSize',14)
    ylabel('Intensity [a.u]', 'FontSize',14)
    hold on
    tick = cellfun(@length, thick_cells);
    a= repmat(OPL_2,[2,1]).*1e6;
    b = a*0;
    b(2,:) = 3e-3;
    plot(a,b, 'LineWidth',2)
     kk = numel(dir(save_path));
    %save([save_path, 'ML (', num2str(kk-2+1), ').mat'],'OPL_2', 'Depth');
end

