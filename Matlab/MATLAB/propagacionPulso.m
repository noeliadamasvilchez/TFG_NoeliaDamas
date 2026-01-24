% =========================================================================
% CÓDIGO BASE
%
% Este script se basa en un código base proporcionado (no desarrollado
% originalmente en el TFG). No se ha modificado la lógica del código;
% únicamente se ha utilizado/ejecutado para generar resultados y figuras.
%
% =========================================================================



% Setup Inicial

tic

clear all;
warning off;

% Parámetros laser pulsado
fRep = 80e6;            % Frecuencia de repetición  [Hz]
FWHM = 0.1;             % Anchura del pulso [pseg]
C = 2;                  % Parámetro Chirp
dutyCycle = 4*FWHM*(1e-12)*fRep/sqrt(2*log(2));
Xi = 4/(dutyCycle*sqrt(pi));
Xi_dB = 10*log10(Xi);

% Physical Parameters
n0 = 1;
c = 3e8;                % Speed of light [m/seg]
wvl = 1550e-9;          % Optical wavelength [m]
k = 2*pi/wvl;           % Optical wavenumber [1/m]
w0 = 2*pi/wvl;          % Angular frequency of light radiation [rad/seg]
f0 = w0/2*pi;           % Frequency of light radiation [Hz]
P0 = 1;                 % Power of the pulse [kW]
fc = 0;                 % Frequency of the carrier [Hz]
wc = 2*pi*fc;           % Angular frequency of the carrier [rad/seg]

% Pulso gaussiano con parámetro chirp C
sigma0 = FWHM/(2*sqrt(2*log(2)));   % Anchura eficaz del pulso [pseg]
F_IN = @(t) gaussianPulse_function(t,sigma0,C);
% Anchura del pulso sobre max de Intensidad/2 en el dominio de la frecuencia
FWHM_frec = (1/2*pi)*sqrt((2*log(2)*(1 + C^2))/sigma0*sigma0);
% Producto anchura temporal x anchura espectral
BWxT_ana = (2/pi)*log(2)*(1 + C^2);
%F_IN = @(t) secantHyperbolic_function(t,sigma0,C);

% Parámetros del canal FSO terrestre
Dz = 3;                 % Distance [km]
h = 500;                % Altura [m]
alpha = 0;              % Atenuación
[Beta0,Beta1,Beta2,Beta3,coefD,longDispersion2,longDispersion3,tauZ] = canalFSO(Dz,sigma0,h,2*pi*c/wvl,c,wvl);

% Pulso transmitido: código independiente de la forma de pulso
tmax = 10*sigma0;                               % Max values of time 
n = 2^12;                                       % Number of points in the numerical mesh of the pulse
tini = linspace(-tmax,tmax,n);                  % Time space
dtini = abs(tini(1)-tini(2));                   % Time space step
dwini = 1/n/dtini*2*pi;                         % Reciprocal space step
omegaini = (-n/2:1:n/2-1)*dwini;                % Frequency space
Carrier = exp(wc*tini*i);
F_in = F_IN(tini);
F0 = sqrt(P0)./sqrt(sum(abs(F_in).^2)*dtini);   % Potencia del pulso P0
F_in = F0*F_in;                                 % Envolvente pulso transmitido
E_in = F_in.*Carrier;                           % Campo eléctrico
fourierF_in = fftshift(fft(F_in));              % TF del pulso en transmisión

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION NLSE %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation parameters
nn = 1e4;                       % Number of propagation steps
tmax = 120*sigma0;              % Max values of time 
t = linspace(-tmax,tmax,n);     % Time space
dt = abs(t(1)-t(2));            % Time space step
dw = 1/n/dt*2*pi;               % Reciprocal space step
omega = (-n/2:1:n/2-1)*dw;      % Frequency space

Carrier_out = exp((wc*t - Beta0*Dz)*i);     % Portadora se propaga a velocidad de fase wc/Beta0
[E_out,F_out,fourierF_out] = simulationNLSE(F0*F_IN(t),Dz,alpha,Beta1,Beta2,Beta3,Carrier_out,omega,n,nn);
P_out = sum(abs(F_out).^2)*dt;
FWHMout_sim = anchuraPulso(F_out,t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Expresión analítica de la propagacion de pulso Gaussiano (asumiento Beta1=0 y por tanto tauZ=0)
[F_out_ana,E_out_ana,sigmaZ0,FWHMout] = propagacionGaussianoAnalitico(t,F0,C,0,Dz,sigma0,Beta2,Carrier_out);

% Representacion pulso gaussiano en transmision
figure(1)
subplot(2,2,1);
plot(tini,real(E_in),'b','LineWidth',1.5);
hold on;
plot(tini,abs(F_in),'r','LineWidth',1.5);
hold off;
title(['Pulso gaussiano en z = 0 m']);
xlabel('Tiempo, t[ps]');ylabel('Amplitud');
legend({'Real(E_{in})','|F_{in}|'},'Location','northwest');
grid on;
ylim([min(ylim) 2.5]);
ax = gca;ax.FontSize = 13;

% Representacion pulso gaussiano tras propagarse una distancia L
subplot(2,2,2);
plot(t,real(E_out),'b','LineWidth',1.5);
hold on;
plot(t,abs(F_out),'r','LineWidth',1.5);
plot(t,real(E_out_ana),'k--','LineWidth',1.5);
plot(t,abs(F_out_ana),'k--','LineWidth',1.5);
hold off;
title(['Pulso gaussiano en z = ',num2str(Dz),' km']);
xlabel('Tiempo, t[ps]');ylabel('Amplitud');
legend({'Real(E_{out})','|F_{out}|','Real(E_{out}) - Ana','|F_{out}| - Ana'},'Location','northwest');
grid on;ax = gca;ax.FontSize = 13;

% Representacion pulso gaussiano en transmision
subplot(2,2,3);
plot(omegaini/2*pi,abs(fourierF_in),'b','LineWidth',1.5);
title(['T. Fourier Gaussian pulse at z = 0 km']);
xlabel('Frequency, f[THz]');ylabel('abs(FFT)');
grid on;ax = gca;ax.FontSize = 13;

% Representacion TF pulso gaussiano en recepción
subplot(2,2,4);
plot(omega/2*pi,abs(fourierF_out(end,:)),'b','LineWidth',1.5);
title(['T. Fourier Gaussian pulse at z = ',num2str(Dz),' km']);
xlabel('Frequency, f[THz]');ylabel('abs(FFT)');
grid on;ax = gca;ax.FontSize = 13;

% Presentación por pantalla
disp(['--------------------------'])
disp(['PARÁMETROS DEL PULSO EN TX'])
disp(['--------------------------'])
disp(['Longitud de onda [nm]: ',num2str(wvl*1e9)])
disp(['Parámetro chirp (C): ',num2str(C)])
disp(['Duración del pulso (FWHM) [pseg]: ',num2str(FWHM)])
disp(['Duración del pulso (FWHM) [pseg]: ',num2str(anchuraPulso(F_in,tini))])
disp(['Anchura eficaz del pulso [pseg]: ',num2str(sigma0)])
disp(['--------------------------'])
disp(['PARÁMETROS DEL ENLACE FSO '])
disp(['--------------------------'])
disp(['Longitud del enlace [km]: ',num2str(Dz)])
disp(['Altura del enlace [m]: ',num2str(h)])
%disp(['Parámetro Beta0 [1/km]: ',num2str(Beta0)])
%disp(['Velocidad de fase [m/seg]: ',num2str(1e3*w0/Beta0)])
disp(['Parámetro Beta1 [ps/km]: ',num2str(Beta1)])
disp(['Retardo de grupo beta1*z [seg]: ',num2str(tauZ*1e-12)])
disp(['Parámetro Beta2 (GVD) [ps^2/km]: ',num2str(Beta2)])
disp(['Parámetro Beta3 (TOD) [ps^3/km]: ',num2str(Beta3)])
disp(['Coeficiente de dispersión D [ps/km*nm]: ',num2str(coefD)])
disp(['Longitud dispersión Beta2 [km]: ',num2str(longDispersion2)])
disp(['Longitud dispersión Beta3 [km]: ',num2str(longDispersion3)])
disp(['--------------------------'])
disp(['PARÁMETROS DEL PULSO EN RX'])
disp(['--------------------------'])
disp(['Duración del pulso (FWHM) [pseg]: ',num2str(FWHMout)])
disp(['Duración del pulso (FWHM) [pseg]: ',num2str(FWHMout_sim)])
disp(['Factor de ensanchamiento: ',num2str(FWHMout_sim/FWHM)])
disp(['Anchura eficaz del pulso [pseg]: ',num2str(sigmaZ0)])

toc