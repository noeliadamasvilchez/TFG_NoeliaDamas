
function [TPS,BWcoh] = TPSVonKarman(Dz,sigma0,Cn2,l0,L0,c);

% Calcula el factor de ensanchamiento temporal por turbulencia (TPS)
% y el ancho de banda de coherencia temporal (BWcoh) usando
% como modelo de turbulencia el modelo von Kármán.
%
% ENTRADAS 
%   Dz     : distancia de propagación [km]
%   sigma0 : anchura temporal [ps]
%   Cn2    : parámetro de estructura del índice de refracción [m^(-2/3)]
%   l0     : escala interna de turbulencia [m]  
%   L0     : escala externa de turbulencia [m]
%   c      : velocidad de la luz [m/s]
%
% SALIDAS
%   TPS   : factor adimensional de ensanchamiento temporal por turbulencia
%   BWcoh : ancho de banda de coherencia temporal [rad/s]
%

 L = Dz * 1e3;                                  % km -> m (longitud del enlace)
 sigma0_s = sigma0 * 1e-12;                     % ps -> s


 k0 = (2*pi) / L0;                               % número de onda espacial asociado a L0 [rad/m]
 a1 = (0.39 * Cn2 * L * (k0^(-5/3))) / (c^2);    % coef. del modelo

 TPS  = 2 * a1 / sigma0_s^2;                     % factor de ensanchamiento temporal (adimensional)    
 BWcoh = 1 / sqrt(a1);                           % ancho de banda de coherencia temporal [rad/s]

end
