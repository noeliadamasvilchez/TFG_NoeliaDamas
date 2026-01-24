
function [a1, TPS, BWcoh, I, kmin, kmax] = A1_OTOPS(out, Dz, c, L0, sigma0)

% Calcula el coeficiente a1 usando el modelo OTOPS como turbulencia y devuelve TPS y BWcoh.
%
% En la formulación clásica de Andrews para atmósfera turbulenta:
%   a1 = (2π² Dz / c²) ∫₀^∞ κ Φ_n(κ) dκ  ≈  0.39 * Cn² * Dz * L₀^(-5/3) / c²
% donde Φ_n(κ) es el espectro de potencia del índice de refracción.
%
% En este trabajo, el espectro atmosférico se sustituye por el espectro
% oceánico OTOPS, por lo que la integral no admite una forma analítica
% (debemos calcular la integral al ser otro espectro, en vez de poder usar
% la aproximación).
%
% Aquí se evalúa numéricamente:
%
%       I = ∫ k * Φ_n(k) dk     entre kmin=1/L₀ y kmax=30/η
%
% El resultado se introduce en:
%       a1 = (2π² * Dz / c²) * I
%
%
%
% ENTRADAS
%   out         : struct de datos obtenidos del script TPS_OTOPS  (incluye
%                 eta y Phi_n_fun)
%   Dz          : longitud del enlace [m]
%   c           : velocidad de la luz [m/s] 
%   L0          : escala externa [m] 
%   sigma0      : ancho temporal inicial del pulso [ps] 
%
% SALIDAS
%   a1          : coeficiente de turbulencia para TPS/BWcoh
%   TPS         : factor de ensanchamiento temporal (adimensional)
%   BWcoh       : ancho de banda coherente [rad/s]
%   I           : integral I = ∫ k * Φ_n(k) dk
%   kmin,kmax   : límites usados [1/m]


% Definición de límites de integración
eta  = out.eta;         % escala de Kolmogorov
kmin = 1./L0;           % límite inferior 
kmax = 30./eta;         % límite superior 

% Realización de la integral
Phi_n = out.Phi_n_fun;               
integrand = @(k) k .* Phi_n(k);        
I = integral(integrand, kmin, kmax, 'ArrayValued', true, 'RelTol', 1e-4);

% Obtención del coeficiente a1
a1 = (2*pi^2 * Dz / c^2) * I;

% Obtención de TPS y BWcoh
sigma0_s = sigma0 * 1e-12;              % ps -> s
TPS      = 2 * a1 / (sigma0_s^2);       % adimensional
BWcoh    = 1 / sqrt(a1);                % rad/s

end
