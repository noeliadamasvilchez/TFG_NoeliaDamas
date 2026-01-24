
function [BF_disp, BF_turb, BF_total, BWcoh] = factorEnsanchado_ocean( ...
    sigma0_ps, Dz_m, C, Beta2_ps2_per_m, Beta3_ps3_per_m, out, L0_m, c)

% Calcula los factores de ensanchamiento (turbulento, dispersivo y total)
% y el ancho de banda coherente (OTOPS) usando unidades coherentes en metros.
%
% ENTRADAS
%   sigma0_ps          : ancho temporal 
%   Dz_m               : longitud de enlace [m]
%   C                  : chirp (adimensional)
%   Beta2_ps2_per_m    : β2 [ps^2/m]
%   Beta3_ps3_per_m    : β3 [ps^3/m]
%   out                : struct de TPS_OTOPS 
%   L0_m               : escala externa [m]
%   c                  : velocidad de la luz [m/s]
%
% SALIDAS
%   BF_disp, BF_turb, BF_total : factores de ensanchamiento (adimensionales)
%   BWcoh                      : ancho de banda coherente [rad/s]

% 1) a1, TPS y BWcoh desde OTOPS 
    
[a1, TPS, BWcoh] = A1_OTOPS(out, Dz_m, c, L0_m, sigma0_ps);

% 2) Ensanchamiento por turbulencia
    
BF_turb = sqrt(1 + TPS);

% 3) Ensanchamiento por dispersión (con β en ps/m y Dz en m)
    
term2 = Beta2_ps2_per_m * Dz_m / (2 * (sigma0_ps^2));
term1 = 1 + C * term2;  
term3 = (Beta3_ps3_per_m * Dz_m) / (4 * (sigma0_ps^3));
BF_disp  = sqrt( term1^2 + term2^2 + 0.5 * (1 + C^2)^2 * term3^2 );

% 4) Ensanchamiento total (turb + disp) 
   
BF_total = sqrt( TPS + term1^2 + term2^2 + 0.5 * (1 + C^2)^2 * term3^2 );

end
