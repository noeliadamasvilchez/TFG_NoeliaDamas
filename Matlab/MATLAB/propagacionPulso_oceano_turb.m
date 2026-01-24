% =========================================================================
%
% Este script evalúa el ensanchamiento temporal de un pulso gaussiano
% ultracorto (con chirp) propagándose en océano, considerando
% simultáneamente dispersión cromática y turbulencia oceánica. 
% Su estructura está basada en el script propagacionPulso_turb.m
%
% =========================================================================



tic;


clear all;
warning off;

%% ===================== 1) PARÁMETROS LASER PULSADO ======================

fRep = 80e6;            % [Hz]
FWHM = 0.2;               % [ps]
C    = 2;               % Chirp
dutyCycle = 4*FWHM*(1e-12)*fRep/sqrt(2*log(2));
Xi = 4/(dutyCycle*sqrt(pi));
Xi_dB = 10*log10(Xi);

%% ===================== 2) PARÁMETROS FÍSICOS (OCÉANO) ===================

n0 = 1;
c = 3e8;                % [m/s]
wvl = 500e-9;           % [m] (OCÉANO)
k = 2*pi/wvl;           % [1/m]
w0 = 2*pi*c/wvl;        % [rad/s]
f0 = w0/2*pi;           % [Hz]
P0 = 1;                 % [kW]
fc = 0;                 % [Hz]
wc = 2*pi*fc;           % [rad/s]


%% ===================== 3) PULSO DE REFERENCIA ===========================

sigma0 = FWHM/(2*sqrt(2*log(2)));   % [ps]
F_IN = @(t) gaussianPulse_function(t, sigma0, C); 

sigma_ref_ps = sigma0;     % [ps]
sigma_ref_fs = sigma0*1e3; % [fs]


%% ===================== 4) PARÁMETROS DEL CANAL OCÉANO ===================

Dz_m = 100;               % [m] distancia
alpha = 0;                % Atenuación

% Coeficientes del agua
coeff.c0 = 1.334200017;
coeff.A  = -1.8279319e9;
coeff.B  = -1.61794879715e28;
coeff.C  = 7.40833353e-34;
coeff.D  = 5.942519e-67;
coeff.E  = 6.2482635e-98;

% Turbulencia (OTOPS): parámetros
args.SP=35; args.T=10; args.Pdbar=20; args.lambda_nm=500;
args.epsilon=1e-8; args.beta0=0.72; args.Xt=1e-10; args.H_grad=-20;
args.lat=0; args.lon=0;
out = TPS_OTOPS(args);

L0_m = 10;   % [m] escala externa


%% ===================== 5) BARRIDO sigma0 ================================

sigma0_fs_full  = linspace(1, 1000, 600);   % [fs]
sigma0_vec_full = sigma0_fs_full * 1e-3;    % [ps]

BF_disp  = zeros(size(sigma0_vec_full));
BF_turb  = zeros(size(sigma0_vec_full));
BF_total = zeros(size(sigma0_vec_full));
BWcoh    = zeros(size(sigma0_vec_full));

Beta1_vec = zeros(size(sigma0_vec_full));
Beta2_vec = zeros(size(sigma0_vec_full));
Beta3_vec = zeros(size(sigma0_vec_full));
coefD_vec = zeros(size(sigma0_vec_full));
tauZ_vec  = zeros(size(sigma0_vec_full));
LD2_vec   = zeros(size(sigma0_vec_full));
LD3_vec   = zeros(size(sigma0_vec_full));

for idx = 1:numel(sigma0_vec_full)

    sigma0 = sigma0_vec_full(idx);  % [ps]

    % Canal OCÉANO (betas en ps/m, longitudes en m)
    [~, Beta1, Beta2, Beta3, coefD, L_D2_m, L_D3_m, tauZ_ps] = ...
        canalOcean(Dz_m, sigma0, coeff, w0, c, wvl);

    % BF combinado OCÉANO (dispersión + turbulencia OTOPS)
    [BF_disp(idx), BF_turb(idx), BF_total(idx), BWcoh(idx)] = ...
        factorEnsanchado_ocean(sigma0, Dz_m, C, Beta2, Beta3, out, L0_m, c);

    Beta1_vec(idx) = Beta1;
    Beta2_vec(idx) = Beta2;
    Beta3_vec(idx) = Beta3;
    coefD_vec(idx) = coefD;
    tauZ_vec(idx)  = tauZ_ps;
    LD2_vec(idx)   = L_D2_m;
    LD3_vec(idx)   = L_D3_m;
end

%% ===================== 6) VALORES PULSO =================================

BFdisp_ref  = interp1(sigma0_fs_full, BF_disp,  sigma_ref_fs, 'linear');
BFturb_ref  = interp1(sigma0_fs_full, BF_turb,  sigma_ref_fs, 'linear');
BFtotal_ref = interp1(sigma0_fs_full, BF_total, sigma_ref_fs, 'linear');
BWcoh_ref   = interp1(sigma0_fs_full, BWcoh,    sigma_ref_fs, 'linear');

Beta1_ref = interp1(sigma0_fs_full, Beta1_vec, sigma_ref_fs, 'linear');
Beta2_ref = interp1(sigma0_fs_full, Beta2_vec, sigma_ref_fs, 'linear');
Beta3_ref = interp1(sigma0_fs_full, Beta3_vec, sigma_ref_fs, 'linear');
coefD_ref = interp1(sigma0_fs_full, coefD_vec, sigma_ref_fs, 'linear');
tauZ_ref  = interp1(sigma0_fs_full, tauZ_vec,  sigma_ref_fs, 'linear');
LD2_ref   = interp1(sigma0_fs_full, LD2_vec,   sigma_ref_fs, 'linear');
LD3_ref   = interp1(sigma0_fs_full, LD3_vec,   sigma_ref_fs, 'linear');

FWHM_disp  = BFdisp_ref  * FWHM;    % [ps]
FWHM_turb  = BFturb_ref  * FWHM;    % [ps]
FWHM_total = BFtotal_ref * FWHM;    % [ps]



%% ===================== 7) GRÁFICA =======================================
figure(1); clf;
set(gcf,'Units','centimeters','Position',[2 2 14 10]);

LW = 1.6;

maskPlot = (sigma0_fs_full >= 0) & (sigma0_fs_full <= 200);

plot(sigma0_fs_full(maskPlot), BF_total(maskPlot), 'b', 'LineWidth', LW); hold on;
plot(sigma0_fs_full(maskPlot), BF_turb(maskPlot),  'r', 'LineWidth', LW);
plot(sigma0_fs_full(maskPlot), BF_disp(maskPlot),  'k--','LineWidth', LW);


if sigma_ref_fs >= 0 && sigma_ref_fs <= 200
    xline(sigma_ref_fs,'k:','LineWidth',1.2);
    plot(sigma_ref_fs, BFtotal_ref,'bo','MarkerFaceColor','b','MarkerSize',5);
    plot(sigma_ref_fs, BFturb_ref ,'ro','MarkerFaceColor','r','MarkerSize',5);
    plot(sigma_ref_fs, BFdisp_ref ,'ko','MarkerFaceColor','k','MarkerSize',5);
end

grid on;
title('Factor de ensanchado (BF)');
xlabel('Anchura eficaz, sigma   \sigma  [fs]');
ylabel('BF');
legend('Total','Turb.','Disp.','Location','northeast');

xlim([0 200]);
ylim([0 1000]);
yticks(0:100:1000);

hold off;


%% ===================== 8) SALIDA POR TERMINAL ===========================

disp('--------------------------')
disp('PARÁMETROS DEL PULSO EN TX')
disp('--------------------------')
disp(['Longitud de onda [nm]: ',num2str(wvl*1e9,'%.2f')])
disp(['Parámetro chirp (C): ',num2str(C)])
disp(['Duración del pulso (FWHM) [ps]: ',num2str(FWHM,'%.2f')])
disp(['Anchura eficaz del pulso sigma0 [ps]: ',num2str(sigma_ref_ps,'%.2f')])
disp(['Anchura eficaz del pulso sigma0 [fs]: ',num2str(sigma_ref_fs,'%.2f')])
disp('--------------------------')

disp('PARÁMETROS DEL ENLACE (OCÉANO)')
disp('--------------------------')
disp(['Distancia [m]: ', num2str(Dz_m,'%.2f')])
disp(['Beta1 [ps/m] (ref): ', num2str(Beta1_ref,'%.2f')])
disp(['Retardo de grupo beta1*z [ps] (ref): ', num2str(tauZ_ref,'%.2f')])
disp(['Beta2 (GVD) [ps^2/m] (ref): ', num2str(Beta2_ref,'%.2e')])
disp(['Beta3 (TOD) [ps^3/m] (ref): ', num2str(Beta3_ref,'%.2e')])
disp(['Coef. dispersión D [ps/(m·nm)] (ref): ', num2str(coefD_ref,'%.2e')])
disp(['Longitud dispersión L_D2 [m] (ref): ', num2str(LD2_ref,'%.2f')])
disp(['Longitud dispersión L_D3 [m] (ref): ', num2str(LD3_ref,'%.2f')])
disp(['L0 [m]: ', num2str(L0_m,'%.2f')])
disp('--------------------------')

disp('PARÁMETROS DEL PULSO EN RX (a partir de BF)')
disp('--------------------------')
disp(['BF_disp   (solo dispersión):   ', num2str(BFdisp_ref,'%.2f')])
disp(['BF_turb   (solo turbulencia):  ', num2str(BFturb_ref,'%.2f')])
disp(['BF_total  (combinado):         ', num2str(BFtotal_ref,'%.2f')])
disp(['FWHM_disp  [ps]: ', num2str(FWHM_disp,'%.2f')])
disp(['FWHM_turb  [ps]: ', num2str(FWHM_turb,'%.2f')])
disp(['FWHM_total [ps]: ', num2str(FWHM_total,'%.2f')])
disp(['BWcoh [rad/s] (ref): ', num2str(BWcoh_ref,'%.2e')])
disp('--------------------------')

toc;
