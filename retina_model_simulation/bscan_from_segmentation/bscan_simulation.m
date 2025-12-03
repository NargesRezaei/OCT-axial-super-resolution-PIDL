% Generating B-Scan of retina
clc
clearvars
close all

%%%%%%%%%%%%%%%%%%%%
%%%       constants
%%%%%%%%%%%%%%%%%%%%
theta0       = 0;
pol          = 'p';
num_smpl     = 150000;
pixel_size   = 2.8e-6;
n_pixel_size = 1e-6;
N0           = 1.33;
z0_R         = 5e-6;
z0_S         = 70e-6;
Ref_scale_f  = .5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% medium parameters
% N & OPL
% Layer names:    BG   RNFL   GCL    IPL   INL    OPL    ONL    ELM    IS     IS/OSJ OS     RPE    BG
N_1            = [N0   1.38,  1.37,  1.374 1.376, 1.39,  1.398, 1.41,  1.425, 1.405, 1.395, 1.38   N0];
DN_cluster     = [0    0.08   0.01   0.015 0.005  0.02   0.003  0.02   0.003  0.2    0.008  0.3    0];
dn_all         = [0    0.03,  0.015, 0.012 0.001  0.015, 0.01,  0.005, 0.01,  0.03,  0.02,  0.08   0];
p_cluster_full = [0    0.55,  0.5,   0.5   0.6,   0.3,   0.85,  0.9,   0.1,   0.95,  0.3,   0.8    0];

cluster_size_range_full = [
    0 0; % BG
    3 5; % NFL
    3 5; % GCL
    2 5; % IPL
    3 5; % INL
    3 5; % OPL
    3 6; % ONL
    1 2; % ELM
    1 3; % IS
    1 2; % ISOSJ
    3 5; % OS
    3 5; % RPE
    0 0; % BG
    ];
env_after_cluster_range_full = [
    0 0; % BG
    0 1; % NFL
    1 2; % GCL
    1 2; % IPL
    1 2; % INL
    1 2; % OPL
    2 3; % ONL
    1 1; % ELM
    1 2; % IS
    1 1; % ISOSJ
    1 2; % OS
    2 3; % RPE
    0 0; % BG
    ];

bright_spot_layers = [1.37, 1.374, 1.39];  % GCL, IPL, OPL
bright_spot_prob_layers = [0.01, 0.01, 0.05];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spectrume parameters

% Source constants
lambda_c  = 840e-9;  %center lambda
Bandwidth = 100e-9;  %spectrum bandwidth
d_lambda  = .1e-9;   %delta lambda
FWHM      = 40e-9;

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

width_f     = 2*(lambda_c - test02(lambda, AmpSpectrum));
dz          = 2 * log(2)/pi* lambda_c^2/width_f; % axial resolution of OCT
l_c         = dz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Define SBW
sbw_Window = linspace(lambda_i, lambda_f, Llambda_sbw);
sbw_convertor = zeros(Llambda_sbw, length(lambda));
for kk = 1 : Llambda_sbw
   sbw_convertor(kk,:) = exp(-((lambda - lambda_sbw(kk)).^2) / (2 * (d_sbw/2)^2));
   sbw_convertor(kk,:) = sbw_convertor(kk,:) / sum(sbw_convertor(kk,:)); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reference arm
ErR = Reference_Mirror(lambda, N0, z0_R, theta0, pol, AmpSpectrum, 'ideal', Ref_scale_f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sample arm
file = 'D:\Narges\DL-MultilayerOCT\Data\Retina_seg3\seg_14040324\_cropno_1250592_3_E copy.mat';
bscn = load(file);
peak = cell(1,length(bscn.nI));
offset = 25;
for jj = 1 : length(bscn.nI)
    q  = find(bscn.nI(:,jj) ~= 0);
    q2 = diff(q)+1;
    peak{jj} = round([q(1)-offset+1; q2]*pixel_size./(n_pixel_size));
end



tic
for jj = 1 : length(bscn.nI)

    pxl_ranges   = peak{jj};
    num_layers   = length(pxl_ranges);
    param_indx   = length(N_1) - num_layers+1;
    N_layer      = [N_1(1) N_1(param_indx:end)];
    dN_clstr     = [DN_cluster(1) DN_cluster(param_indx:end)];
    p_cluster_s  = [p_cluster_full(1) p_cluster_full(param_indx:end)];
    cluster_size = [cluster_size_range_full(1,:); cluster_size_range_full(param_indx:end, :)];
    dn_all_s     = [dn_all(1) dn_all(param_indx:end)];
    env_s        = [env_after_cluster_range_full(1,:) ;env_after_cluster_range_full(param_indx:end,:)];

    thick_cells  = cell(1, num_layers);
    refidx_cells = cell(1, num_layers);
    OPL          = zeros(1, num_layers);

    thick_cells{1}  = n_pixel_size * ones(1, pxl_ranges(1));
    refidx_cells{1} = N_layer(1)   * ones(1, pxl_ranges(1));
    OPL(1)          = sum(thick_cells{1} .* refidx_cells{1});

    for i = 2:num_layers
        len                = pxl_ranges(i);
        base_value         = N_layer(i);
        base_cluster_value = dN_clstr(i);
        dn                 = dn_all_s(i);
        p_cluster          = p_cluster_s(i);
        cluster_range      = cluster_size(i, :);
        env_range          = env_s(i, :);


        ref_index   = zeros(1, len);
        trans_range = min(3, len);
        prev_base   = N_layer(i - 1);

        % Skip transition for first and last layers
        if i == 1 || i == 2 || i == length(N_layer)
            trans_range = 0; % No transition for first or last layer
        else
            % Smooth transition for middle layers
            for k = 1:trans_range
                alpha = k / trans_range;
                ref_index(k) = (1 - alpha) * prev_base + alpha * base_value + (2 * rand - 1) * dn / 4;
            end
        end

        % Start cursor after transition
        cursor = trans_range + 1;

        while cursor <= len
            if rand < p_cluster
                clen            = min(randi(cluster_range), len - cursor + 1);
                prob_index      = find(abs(base_value - bright_spot_layers) < 1e-4, 1);
                cluster_val     = base_value + (2 * rand - 1) * dn / 4 + base_cluster_value;

                if ~isempty(prob_index) && rand < bright_spot_prob_layers(prob_index)
                    cluster_val     = cluster_val + 0.1;
                end

                % % Add random phase shift for ONL
                % if base_value == 1.385 || base_value == 1.405 || base_value == 1.38
                %     cluster_val = cluster_val + (rand - 0.5) * 0.01; % Random phase
                % end

                ref_index(cursor:cursor+clen-1) = cluster_val ;
                cursor                          = cursor + clen;
                if cursor > len, break; end
                % Add one environment pixel after the cluster
     
                % Normal environment pixel
                env_val = base_value;
                ref_index(cursor) = env_val;
                cursor = cursor + 1;
                if cursor > len, break; end
            else

                n_env = randi(env_range);
                for e = 1:n_env
                    if cursor > len, break; end
                    env_val = base_value ;
                    ref_index(cursor) = env_val;
                    cursor = cursor + 1;
                end
            end

        end

        % if abs(mean(ref_index) - base_value) > 0.001
        %     ref_index = ref_index + (base_value - mean(ref_index));
        % end

        thick_cells{i}  = n_pixel_size * ones(1, len);
        refidx_cells{i} = ref_index;
        OPL(i)          = sum(thick_cells{i} .* ref_index);
    end

    % Concatenate all thickness and refractive index values
    LD = cell2mat(thick_cells);         % total layer thickness array
    LN = cell2mat(refidx_cells);       % total refractive index array

    % Append final refractive index value (bottom layer)
    LN = [LN, N_1(end)];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spectral Interferogram
    [Er, ~, ~, ~, ~] = General_Multilayer_V11(lambda, [N0, LN],  [z0_S  LD], theta0, pol, AmpSpectrum/sqrt(2));
    E_sum = Er +  ErR;
    I_OCT2 = 1/2 * (E_sum .* conj(E_sum));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spectrometer

    I_OCT3 = sbw_convertor * I_OCT2';

    ref2 = 1/2 * (Er .* conj(Er) + ErR .* conj(ErR));
    ref3 = sbw_convertor * ref2';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spectral Interferogram Analyse
    [Xaxis, Depth] = OCT_Analyse(I_OCT3, ref3, lambda_sbw , LD);
    Xaxis = round(Xaxis', 7);
    Depth = Depth(1:200)';
    Xaxis = Xaxis(1:200);
    
    OPL_2 = cumsum(OPL);


    B_Depth(:,jj) = Depth;
end
toc

figure;
bscan_db = 10 * log10(B_Depth + 1e-5);
imagesc(1:612*4, Xaxis*1e6, bscan_db)
colormap gray
axis image



