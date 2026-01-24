function FWHM = anchuraPulso(F,t)

% Anchura temporal
FWHM = 0;
Intensidad = abs(F).^2;
MAX = max(Intensidad);
%error = 5e-3;
%pos_MAX_mitad = find(abs(Intensidad - MAX*0.5) < error);

% Bucle que busca el FWHM
error = 0.01;
long_MAX_mitad = 10;
while long_MAX_mitad > 2
    pos_MAX_mitad = find(abs(Intensidad - MAX*0.5) < error);
    long_MAX_mitad = length(pos_MAX_mitad);
    error = error - 0.00001;
end

if length(pos_MAX_mitad) == 2
    FWHM = 2*t(pos_MAX_mitad(2));
end

% Anchura espectral
% 
% DEE = abs(fft(F(tiempo))).^2;
% 
% MAX = max(DEE);
% pos_MAX_mitad = find(abs(DEE - MAX*0.5) < 0.1);
% 
% if length(pos_MAX_mitad) == 2
%     FWHM_f = 2*(pos_MAX_mitad(2))/tiempo;
% else
%     FWHM_f = 0;
% end