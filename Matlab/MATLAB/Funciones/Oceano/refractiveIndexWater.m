
function n = refractiveIndexWater(wvl, coeff)

% Índice de refracción del agua salada (óceano)
%
% ENTRADAS
%   wvl  : longitud de onda [m]
%   coeff: struct con campos A, B, C, D, E, c0 (evaluados a T, p, S)
%          n(ω) = c0 + A/ω + B/ω^2 + C*ω^2 + D*ω^4 + E*ω^6
% 
% SALIDAS
%   n    : índice de refracción
%

c = 3e8;                    
omega = 2*pi*c ./ wvl;          
A = coeff.A; B = coeff.B; C = coeff.C; D = coeff.D; E = coeff.E; c0 = coeff.c0;
n = c0 + A./omega + B./(omega.^2) + C.*(omega.^2) + D.*(omega.^4) + E.*(omega.^6);

end
