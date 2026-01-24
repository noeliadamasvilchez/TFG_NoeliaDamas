

function n = refractiveIndex(wvl,h)

% h altura en metros
% wvl longitud de ondas en m

% Constantes
a = 77.6e-6;
b = 7.52e-15;

% Temperatrura y presión
Th = 288.19 - h*6.49e-3;                % Kelvin
Ph = ((44.41 - h*1e-3)^5.256)*2.23e-6;  % mbar

n = 1 + a*Ph*(1/Th)*(1 + b*(wvl.^(-2))); % Índice refacción Andrews