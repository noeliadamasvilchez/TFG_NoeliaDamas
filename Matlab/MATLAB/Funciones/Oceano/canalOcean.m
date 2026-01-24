
function [Beta0,Beta1,Beta2,Beta3,coefD,longDispersion2,longDispersion3,tauZ] = ...
         canalOcean(Dz, sigma0, coeff, w, c, wvl)

% Versión oceánica adaptada a unidades en metros y betas en ps/m del canal
%
% ENTRADAS
%   Dz     [m]        longitud del enlace
%   sigma0 [ps]       ancho temporal
%   coeff  (struct)   coeficientes del modelo de índice del agua
%   w      [rad/s]    frecuencia angular central
%   c      [m/s]      velocidad de la luz
%   wvl    [m]        longitud de onda
%
% SALIDAS
%   Beta0  [1/m]
%   Beta1  [ps/m]
%   Beta2  [ps^2/m]
%   Beta3  [ps^3/m]
%   coefD  [ps/(m·nm)]
%   longDispersion2 [m]
%   longDispersion3 [m]
%   tauZ   [ps]

% 1) Betas en SI desde el modelo de índice del agua 

[Beta0_SI, Beta1_SI, Beta2_SI, Beta3_SI] = taylorBeta_ocean(coeff, w, c);

% 2) Coeficiente de dispersión D en SI y luego a ps/(m·nm) 
   
 coefD_SI = -2*pi*c*Beta2_SI/(wvl*wvl);   % [s/m^2]
 coefD    = coefD_SI * 1e3;               % [ps/(m·nm)]  

 % 3) Conversión de betas a ps/m, ps^2/m, ps^3/m 

 Beta0 = Beta0_SI;           % [1/m]
 Beta1 = Beta1_SI * 1e12;    % [s/m]  -> [ps/m]
 Beta2 = Beta2_SI * 1e24;    % [s^2/m]-> [ps^2/m]
 Beta3 = Beta3_SI * 1e36;    % [s^3/m]-> [ps^3/m]


 % 4) Longitudes de dispersión (en metros)
 % Sustituimos sigma0 por su equivalente en FWHM directamente 
 % en las fórmulas -> FWHM = sigma0 * (2 * sqrt(2 * log(2)))

 % Fórmula teoría LD2: FWHM^2 / (4 * ln(2) * |Beta2|)
 longDispersion2 = ( (sigma0 * 2 * sqrt(2 * log(2)))^2 ) / (4 * log(2) * abs(Beta2)); 

 if Beta3 == 0
     longDispersion3 = Inf;
 else
     % Fórmula teoría LD3: FWHM^3 / (8 * (ln(2))^2 * |Beta3|)
     longDispersion3 = ( (sigma0 * 2 * sqrt(2 * log(2)))^3 ) / (8 * (log(2))^2 * abs(Beta3));
 end
 

 % 5) Retardo de grupo 
 tauZ = Beta1 * Dz;   % [ps]  
 Beta0 = 0;
 Beta1 = 0;

end
