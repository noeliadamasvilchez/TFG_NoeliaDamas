
function F = gaussianPulse_function(t,sigma0,C)

% Pulso gaussiano con parámetro chirp C
% t variable temporal
% F0 amplitud
% sigma0 = 1;       % Anchura eficaz del pulso
% C = 0;            % Parámetro Chirp

a = 1/(4*sigma0*sigma0);
b = C/(4*sigma0*sigma0);

% Expresión del pulso en el dominio del tiempo
F = exp(-(a - b*i)*t.*t);