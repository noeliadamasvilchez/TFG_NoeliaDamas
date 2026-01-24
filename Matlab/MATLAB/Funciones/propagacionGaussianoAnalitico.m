

function [F_out,E_out,sigmaZ0,FWHM] = propagacionGaussianoAnalitico(t,F0,C,tauZ,Dz,sigma0,Beta2,Carrier_out)

% Expresión analítica de la propagacion de pulso Gaussiano
sigmaZ0 = sqrt((sigma0 + 0.5*C*Beta2*Dz/sigma0)^2 + (0.5*Beta2*Dz/sigma0)^2);
FWHM = sigmaZ0*(2*sqrt(2*log(2)));
sigma0Chirp = sqrt((sigma0^2)/(1 - C*i));
sigmaZ = sqrt(sigma0Chirp^2 + 0.5*Beta2*Dz*i);
F_out = F0*(sigma0Chirp/sigmaZ).*exp(-((t - tauZ).^2)/(4*sigmaZ*sigmaZ));
E_out = F_out.*Carrier_out;