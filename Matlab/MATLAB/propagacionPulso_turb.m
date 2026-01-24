% =========================================================================
%
% Este script corresponde a un código utilizado para simular la
% propagación de un pulso gaussiano con chirp a través de un enlace FSO,
% considerando los efectos de dispersión cromática y turbulencia atmosférica.
% Está basado en el script propagacionPulso.m
%
% =========================================================================

tic


clear all;
warning off;

%% ===================== 1) PARÁMETROS LASER PULSADO ======================

fRep = 80e6;            % [Hz]
FWHM = 2;               % [ps]
C = 2;                  % Chirp
dutyCycle = 4*FWHM*(1e-12)*fRep/sqrt(2*log(2));
Xi = 4/(dutyCycle*sqrt(pi));
Xi_dB = 10*log10(Xi);


%% ===================== 2) PARÁMETROS FÍSICOS ============================

n0 = 1;
c = 3e8;                % [m/s]
wvl = 1550e-9;          % [m]
k = 2*pi/wvl;           % [1/m]
w0 = 2*pi*c/wvl;        % [rad/s]
f0 = w0/2*pi;           % [Hz]
P0 = 1;                 % [kW] 
fc = 0;                 % [Hz]
wc = 2*pi*fc;           % [rad/s]


%% ===================== 3) PULSO DE REFERENCIA ===========================

sigma0 = FWHM/(2*sqrt(2*log(2)));   % [ps] 
F_IN = @(t) gaussianPulse_function(t,sigma0,C); 
sigma_ref_ps = sigma0;          % [ps]
sigma_ref_fs = sigma0*1e3;      % [fs]


%% ===================== 4) PARÁMETROS DEL CANAL FSO ======================

Dz = 3;                 % [km]
h  = 500;               % [m]
alpha = 0;              % Atenuación

% Turbulencia (se cambia el parámetro Cn2 para simular régimen débil/moderado/fuerte de turbulencia atmosférica)
Cn2 = 1e-12;            % [m^(-2/3)]
l0  = (1e-9 * h)^(1/3); % [m]
L0  = sqrt(4*h);        % [m]



%% ===================== 5) BARRIDO sigma0 ================================

sigma0_fs  = linspace(1, 200, 400);     % [fs]
sigma0_vec = sigma0_fs * 1e-3;          % [ps]

BF_disp  = zeros(size(sigma0_vec));
BF_turb  = zeros(size(sigma0_vec));
BF_total = zeros(size(sigma0_vec));
BWcoh    = zeros(size(sigma0_vec));

Beta1_vec = zeros(size(sigma0_vec));
Beta2_vec = zeros(size(sigma0_vec));
Beta3_vec = zeros(size(sigma0_vec));
coefD_vec = zeros(size(sigma0_vec));
tauZ_vec  = zeros(size(sigma0_vec));
LD2_vec   = zeros(size(sigma0_vec));
LD3_vec   = zeros(size(sigma0_vec));

for idx = 1:numel(sigma0_vec)

    sigma0 = sigma0_vec(idx);  % [ps]

    [~,Beta1,Beta2,Beta3,coefD,longDispersion2,longDispersion3,tauZ] = ...
        canalFSO(Dz,sigma0,h,2*pi*c/wvl,c,wvl);

    [BF_disp(idx), BF_turb(idx), BF_total(idx), BWcoh(idx)] = ...
        factorEnsanchado(sigma0,Dz,C,Beta2,Beta3,Cn2,l0,L0,c);

    Beta1_vec(idx) = Beta1;
    Beta2_vec(idx) = Beta2;
    Beta3_vec(idx) = Beta3;
    coefD_vec(idx) = coefD;
    tauZ_vec(idx)  = tauZ;
    LD2_vec(idx)   = longDispersion2;
    LD3_vec(idx)   = longDispersion3;
end


%% ===================== 6) VALORES PULSO =================================
BFdisp_ref  = interp1(sigma0_fs, BF_disp,  sigma_ref_fs, 'linear', 'extrap');
BFturb_ref  = interp1(sigma0_fs, BF_turb,  sigma_ref_fs, 'linear', 'extrap');
BFtotal_ref = interp1(sigma0_fs, BF_total, sigma_ref_fs, 'linear', 'extrap');
BWcoh_ref   = interp1(sigma0_fs, BWcoh,    sigma_ref_fs, 'linear', 'extrap');

Beta1_ref = interp1(sigma0_fs, Beta1_vec, sigma_ref_fs, 'linear', 'extrap');
Beta2_ref = interp1(sigma0_fs, Beta2_vec, sigma_ref_fs, 'linear', 'extrap');
Beta3_ref = interp1(sigma0_fs, Beta3_vec, sigma_ref_fs, 'linear', 'extrap');
coefD_ref = interp1(sigma0_fs, coefD_vec, sigma_ref_fs, 'linear', 'extrap');
tauZ_ref  = interp1(sigma0_fs, tauZ_vec,  sigma_ref_fs, 'linear', 'extrap');
LD2_ref   = interp1(sigma0_fs, LD2_vec,   sigma_ref_fs, 'linear', 'extrap');
LD3_ref   = interp1(sigma0_fs, LD3_vec,   sigma_ref_fs, 'linear', 'extrap');

FWHM_disp  = BFdisp_ref  * FWHM;    % [ps]
FWHM_turb  = BFturb_ref  * FWHM;    % [ps]
FWHM_total = BFtotal_ref * FWHM;    % [ps]


%% ===================== 7) GRÁFICA =======================================

figure(1);
set(gcf,'Units','centimeters','Position',[2 2 14 10]); 

LW = 1.6;

plot(sigma0_fs, BF_total, 'b', 'LineWidth', LW); hold on;
plot(sigma0_fs, BF_turb,  'r', 'LineWidth', LW);
plot(sigma0_fs, BF_disp,  'k--','LineWidth', LW);

if sigma_ref_fs >= 0 && sigma_ref_fs <= 200
    xline(sigma_ref_fs,'k:','LineWidth',1.2);
    plot(sigma_ref_fs, BFtotal_ref,'bo','MarkerFaceColor','b','MarkerSize',5);
    plot(sigma_ref_fs, BFturb_ref ,'ro','MarkerFaceColor','r','MarkerSize',5);
    plot(sigma_ref_fs, BFdisp_ref ,'ko','MarkerFaceColor','k','MarkerSize',5);
end

hold off;
grid on;
title('Factor de ensanchado (BF)');
xlabel('Anchura eficaz, sigma   \sigma  [fs]');
ylabel('BF');
legend('Total','Turb.','Disp.','Location','northeast');

xlim([0 200]);
ylim([0 1000]);        
yticks(0:100:1000);


%% ===================== 8) SALIDA POR TERMINAL ===========================

disp('--------------------------')
disp('PARÁMETROS DEL PULSO EN TX')
disp('--------------------------')
disp(['Longitud de onda [nm]: ', num2str(wvl*1e9, '%.2f')])
disp(['Parámetro chirp (C): ', num2str(C)])
disp(['Duración del pulso (FWHM) [ps]: ', num2str(FWHM, '%.2f')])
disp(['Anchura eficaz del pulso sigma0 [ps]: ', num2str(sigma_ref_ps, '%.2f')])
disp(['Anchura eficaz del pulso sigma0 [fs]: ', num2str(sigma_ref_fs, '%.2f')])
disp('--------------------------')

disp('PARÁMETROS DEL ENLACE FSO ')
disp('--------------------------')
disp(['Longitud del enlace [km]: ', num2str(Dz, '%.2f')])
disp(['Altura del enlace [m]: ', num2str(h, '%.0f')])
disp(['Cn2 [m^(-2/3)]: ', num2str(Cn2,'%.0e')])
disp(['L0  [m]: ', num2str(L0, '%.2f')])
disp(['l0  [m]: ', num2str(l0, '%.2f')])
disp(['Parámetro Beta1 [ps/km] (ref): ', num2str(Beta1_ref, '%.2f')])
disp(['Retardo de grupo beta1*z [s] (ref): ', num2str(tauZ_ref*1e-12, '%.2e')])
disp(['Parámetro Beta2 (GVD) [ps^2/km] (ref): ', num2str(Beta2_ref, '%.2f')])
disp(['Parámetro Beta3 (TOD) [ps^3/km] (ref): ', num2str(Beta3_ref, '%.2f')])
disp(['Coeficiente de dispersión D [ps/(km·nm)] (ref): ', num2str(coefD_ref, '%.2f')])
disp(['Longitud dispersión Beta2 [km] (ref): ', num2str(LD2_ref, '%.2f')])
disp(['Longitud dispersión Beta3 [km] (ref): ', num2str(LD3_ref, '%.2f')])
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

toc
