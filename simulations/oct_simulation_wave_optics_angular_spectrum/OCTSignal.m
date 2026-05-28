%%% OCT Signal Simulator


clc
clearvars
close all

%%%%%%%%%%%%%%%%%%%%
%%%       constants
%%%%%%%%%%%%%%%%%%%%
tic
% lambda constants
lambda_c  = 850e-9;  %center lambda
Bandwidth = 240e-9;  %spectrum bandwidth
d_lambda  = .1e-9;    %delta lambda

% Incident angle in rad.
lambda0   = 550e-9; % Lambda for thickness of layers.

% medium parameters
z0 = 10; %distance from source
N0 = 1; %refractive index of our medium
%Ni = [1.38     , 2.02     , 1.8      ,1.5] ; %refractive indexes of our sample
%Di = [lambda0/4, lambda0/4, lambda0/4]; % lengths of  our sample
%for g = 1 : 100

dz = log(2)/pi* lambda_c^2/Bandwidth; % axial resolution of OCT

%Ni = N0 + 2*rand(1,4);
%Di =  10*dz + 2*dz*rand(1,3);
Num_reflections = 4; %how many reflection we want to calculate

Ni = [1 1.5 1.2 1];
Di =  [300e-6 250e-6 150e-6];

% space constants
Nx = 1; %number of pixels
Ny = 1;

x0_min = -5e-3;
x0_max =  5e-3;
y0_min = -5e-3;
y0_max =  5e-3;

wg = 200e-3;
%%%%%%%%%%%%%%%%%%%%
%       calculations
%%%%%%%%%%%%%%%%%%%%

%Thickness
%Di = Di./Ni(1:end-1);

% lambda
lambda_i = lambda_c - Bandwidth/2;
lambda_f = lambda_c + Bandwidth/2;
Llambda = Bandwidth /d_lambda+1;
lambda = linspace( lambda_i, lambda_f, Llambda - mod(Llambda, 2) );
AmpSpectrum =  gausswin( length(lambda), 3); %amplitude
%AmpSpectrum =  lambda*0+1; %amplitude

% space
x0 = linspace(x0_min,x0_max,Nx);
y0 = linspace(y0_min,y0_max,Ny);
[x,y] = meshgrid(x0,y0);
% dxo = x0(2) - x0(1);
% dyo = y0(2) - y0(1);
r = sqrt(x.^2+y.^2);
A = exp(-(r/wg).^2);
N.Ni = Ni;
N.N0 = N0;
Ini.rfl = Num_reflections;

%%%%%%%%%%%%%%%%%%%%
%       Propagation
%%%%%%%%%%%%%%%%%%%%

% 1st arm (sample arm)
U_S = zeros(Nx, Ny, length(lambda));
OPL = N0*sqrt(x.^2+y.^2+z0^2);

for ii = 1 : length(lambda)
    Ini.lambda = lambda(ii);
    U.Ui = AmpSpectrum(ii).*A.*exp(1i.*OPL*2*pi/Ini.lambda); %initial wave for each wavelengths
    U_S(:, :, ii) = OCTsimulator(U, N, Di, Ini);
    %     imagesc(abs(U_S(:, :, ii)).^2);
    %     title(num2str(ii))
    %     pause(.1)

end

% 2nd arm (refrence arm)
for ii = 1:length(lambda)
    U_r(:,:,ii) =0.5.*A.*AmpSpectrum(ii).*exp(1i.*(2*pi/lambda(ii))*(x*0+1)*(N0*z0));
    %    imagesc(abs(U_r(:, :, ii)).^2);
    %    title(num2str(ii))
    %    pause(.1)
end

% Interference
U_sum = U_S + U_r;
I_sum = U_sum .* conj(U_sum);
I_OCT =  smooth(squeeze(I_sum(ceil(end/2),ceil(end/2),:)));

figure
plot(lambda, I_OCT)

I_1 = U_r .* conj(U_r);
I_1 = squeeze(I_1(ceil(end/2),ceil(end/2),:));

I_2 = U_S .* conj(U_S);
I_2 = squeeze(I_2(ceil(end/2),ceil(end/2),:));

ref  = I_1 + I_2;

%Analyse and plot
figure
[ Xaxis, Depth ] = OCT_Analyse(I_OCT, ref, lambda, Di);
Depth = Depth/max(Depth);
plot(Xaxis, Depth);
d = findpeaks(Depth);
aa = d(d > (1/15));
for ii = 1 : numel(aa)
    ad(ii) = Xaxis(Depth == aa(ii));
end
ad = [0 ad];
disp(diff(ad));
disp(Di.*Ni(1:end-1));

%hold on
%plot(lambda, I_2,'color', 'g');
toc
%plot(Xaxis, Depth)
%figure, plot(I_OCT)
