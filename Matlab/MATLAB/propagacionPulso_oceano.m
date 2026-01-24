% =========================================================================
% 
% Este script simula la propagación de un pulso gaussiano ultracorto (con
% chirp) a través de un enlace óptico en medio oceánico considerando
% únicamente los efectos de dispersión cromática. Su estructura está basada
% en el script propagacionPulso.m pero adaptado al nuevo medio oceánico.

% =========================================================================

tic; 

clear all; 
warning off;

% Parametros laser pulsado (igual que en atm)
fRep = 80e6;            % Frecuencia de repetición  [Hz]
FWHM = 0.2;             % Anchura del pulso [pseg]
C = 2;                  % Parámetro Chirp
dutyCycle = 4*FWHM*(1e-12)*fRep/sqrt(2*log(2));
Xi = 4/(dutyCycle*sqrt(pi));
Xi_dB = 10*log10(Xi);


% Physical Parameters
n0 = 1;
c = 3e8;                % Speed of light [m/seg]
wvl = 500e-9;           % Optical wavelength [m] CAMBIA RESPECTO ATM
k = 2*pi/wvl;           % Optical wavenumber [1/m]
w0 = 2*pi*c/wvl;        % Angular frequency of light radiation [rad/seg]
f0 = w0/2*pi;           % Frequency of light radiation [Hz]
P0 = 1;                 % Power of the pulse [kW]
fc = 0;                 % Frequency of the carrier [Hz]
wc = 2*pi*fc;           % Angular frequency of the carrier [rad/seg]


% Pulso gaussiano con parámetro chirp C
sigma0 = FWHM/(2*sqrt(2*log(2)));   % Anchura eficaz del pulso [pseg]
F_IN = @(t) gaussianPulse_function(t, sigma0, C);
% Anchura del pulso sobre max de Intensidad/2 en el dominio de la frecuencia
FWHM_frec = (1/2*pi)*sqrt((2*log(2)*(1 + C^2))/sigma0*sigma0);
% Producto anchura temporal x anchura espectral
BWxT_ana = (2/pi)*log(2)*(1 + C^2);


% Parámetros del canal oceánico (cambia respecto atm)
Dz_m = 100;                  % Distancia [m] 
alpha = 0;                   % Atenuación


% Coeficientes del agua (T, p, S evaluados)
coeff.c0 = 1.334200017;
coeff.A  = -1.8279319e9;
coeff.B  = -1.61794879715e28;
coeff.C  = 7.40833353e-34;
coeff.D  = 5.942519e-67;
coeff.E  = 6.2482635e-98;

[Beta0, Beta1, Beta2, Beta3, coefD, L_D2_m, L_D3_m, tauZ_ps] = ...
    canalOcean(Dz_m, sigma0, coeff, w0, c, wvl);


% Adaptación de malla temporal
[t, dt, omega, sigmaZ_est] = auto_time_grid_m(sigma0, Beta2, Dz_m, C, 2^18, 2^19, 8);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION NLSE %%%%%%%%%%%%%%%%%%%%%%%%%%%
Ftx = F_IN(t);                                      % Evaluación del pulso de entrada sobre malla temporal adaptada
F0  = sqrt(P0)/sqrt(sum(abs(Ftx).^2)*dt);           % Factor de normalización
Ftx = F0*Ftx;                                       % Normalización pulso de entrada
Carrier     = exp(1i*wc*t);                         % Portadora en tx
Carrier_out = exp(1i*(wc*t - Beta0*Dz_m));          % Portadora en rx 
nn = 500;                                           % Pasos de propagación (500 para 0.1 de FWHM, sino 1e3)

[E_out, F_out, fourierF_out] = simulationNLSE( ...
    Ftx, Dz_m, alpha, Beta1, Beta2, Beta3, Carrier_out, omega, numel(t), nn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Expresión analítica de la propagacion de pulso Gaussiano
[F_out_ana, E_out_ana, sigmaZ0_ps, FWHMout_ana] = ...
    propagacionGaussianoAnalitico(t, F0, C, 0, Dz_m, sigma0, Beta2, Carrier_out);


% Representacion pulso gaussiano 
figure(1); 
clf;

% Gráfica pulso en tx y rx en tiempo
subplot(2,1,1);
plot(t, abs(Ftx)/max(abs(Ftx)), 'b-', 'LineWidth', 1.2);                % Pulso en tx en tiempo 
hold on;
plot(t, abs(F_out)/max(abs(F_out)), 'r-', 'LineWidth', 1.5);            % Pulso en rx en tiempo
plot(t, abs(F_out_ana)/max(abs(F_out_ana)), 'k:', 'LineWidth', 2);      % Pulso en rx analítico
grid on;
xlabel('Tiempo, t[ps]'); ylabel('Amplitud normalizada');
title(['Pulso gaussiano en z = ', num2str(Dz_m), ' m']);
legend('|F_{in}|', '|F_{out}|', '|F_{out}| - Ana', 'Location', 'northwest');
xlim([-8*sigmaZ0_ps, 8*sigmaZ0_ps]); 

% Gráfica pulso en tx y rx en frecuencia
subplot(2,1,2);
plot(omega/(2*pi), abs(fftshift(fft(Ftx))), 'b');                       % Pulso en tx en freq
hold on;
plot(omega/(2*pi), abs(fourierF_out(end,:)), 'r--');                    % Pulso en rx en freq
grid on;
xlabel('Frecuencia [THz]'); ylabel('Magnitud');
title('Transformada de Fourier');
legend('Espectro TX', 'Espectro RX');


FWHMtx_sim  = anchuraPulso(Ftx, t);
FWHMout_sim = anchuraPulso(F_out, t);


% Presentación por terminal
disp('--------------------------')
disp('PARÁMETROS DEL PULSO EN TX')
disp('--------------------------')
disp(['FWHM inicial [ps]: ', num2str(FWHM)])
disp(['FWHM medido TX [ps]: ', num2str(FWHMtx_sim)])
disp('--------------------------')
disp('PARÁMETROS DEL ENLACE (OCÉANO)')
disp('--------------------------')
disp(['Distancia [m]: ', num2str(Dz_m)])
disp(['Beta2 [ps^2/m]: ', num2str(Beta2)])
disp(['Beta3 [ps^3/m]: ', num2str(Beta3)])
disp(['Longitud Dispersión L_D2 [m]: ', num2str(L_D2_m)])
disp(['Longitud Dispersión L_D3 [m]: ', num2str(L_D3_m)])
disp('--------------------------')
disp('RESULTADOS DE ENSANCHAMIENTO')
disp('--------------------------')
disp(['FWHM Analítico RX [ps]: ', num2str(FWHMout_ana)])
disp(['FWHM Simulado RX  [ps]: ', num2str(FWHMout_sim)])
disp(['FACTOR DE ENSANCHAMIENTO: ', num2str(FWHMout_sim/FWHMtx_sim)])

toc

% Función auxiliar: crea una malla temporal/espectral adecuada para simular el pulso
% teniendo en cuenta el ensanchamiento esperado por dispersión y el chirp (C), siendo
% mucho más acusado que en el caso de atmósfera
function [t, dt, omega, sigmaZ_est] = auto_time_grid_m(sigma0_ps, Beta2_ps_per_m, Dz_m, C, n_min, n_max, margin_sigmas)
    if nargin<6, n_max = 2^19; end
    if nargin<5, n_min = 2^18; end
    if nargin<7, margin_sigmas = 6 ; end        % margin_sigmas = 6 para 0.1 de FWHM

    term = (Beta2_ps_per_m * Dz_m) / (2*sigma0_ps);        
    sigmaZ_est = sqrt( (sigma0_ps + C*term).^2 + term.^2 ); 
    
    tmax = margin_sigmas * sigmaZ_est;        
  
    dt_target = min(sigma0_ps/10, sigmaZ_est/256);               
    
    n = 2^nextpow2( ceil(2*tmax/dt_target) );
    n = min(max(n, n_min), n_max);
  
    
    t = linspace(-tmax, tmax, n);             
    dt = t(2)-t(1);                           
    dw = 2*pi/(n*dt);                         
    omega = (-n/2:n/2-1)*dw;                  
end

   