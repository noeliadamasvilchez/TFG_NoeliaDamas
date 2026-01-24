
function [Beta0,Beta1,Beta2,Beta3] = taylorBeta_ocean(coeff, w0, c)

%
% ENTRADAS
%   coeff : struct con A,B,C,D,E,c0 para n(ω) vía refractiveIndexWater
%   w0    : frecuencia angular central [rad/s]
%   c     : velocidad de la luz [m/s]
%


% Se ha definido w como variable simbólica real (frecuencia angular) para
% obtener β(ω) y sus derivadas dβ/dω, d²β/dω², d³β/dω³ de forma exacta.
% n(ω) se obtiene a partir de n(λ) usando λ = 2πc/ω, de modo que:
%   β(ω) = (ω/c) · n(λ(ω))
syms w real
n_of_w = @(ww) refractiveIndexWater( 2*pi*c./ww, coeff );  
betaSym = (w/c) * n_of_w(w);                               

d1 = diff(betaSym, w, 1);   
d2 = diff(betaSym, w, 2);   
d3 = diff(betaSym, w, 3);   

Beta0 = double( subs(betaSym, w, w0) );   % [1/m]
Beta1 = double( subs(d1,      w, w0) );   % [s/m]
Beta2 = double( subs(d2,      w, w0) );   % [s^2/m]
Beta3 = double( subs(d3,      w, w0) );   % [s^3/m]

end
