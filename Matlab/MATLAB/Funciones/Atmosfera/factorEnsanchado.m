


function [BF_disp,BF_turb,BF_total,BWcoh] = factorEnsanchado(sigma0,Dz,C,Beta2,Beta3,Cn2,l0,L0,c)

[TPS,BWcoh] = TPSVonKarman(Dz,sigma0,Cn2,l0,L0,c);

BF_turb = sqrt(1+TPS); % adimensional, es el factor de ensanchamiento
BF_disp = sqrt((1 + Beta2*Dz*C/(2*(sigma0^2)))^2 + (Beta2*Dz/(2*(sigma0^2)))^2 + 0.5*((1 + C^2)^2)*((Beta3*Dz/(4*(sigma0^3)))^2)); % también adimensional

% De momento solamente es válida cuando C=0 y Beta3=0
BF_total = sqrt(TPS + (1 + Beta2*Dz*C/(2*(sigma0^2)))^2 + (Beta2*Dz/(2*(sigma0^2)))^2 + 0.5*((1 + C^2)^2)*((Beta3*Dz/(4*(sigma0^3)))^2)); % también adimensional
