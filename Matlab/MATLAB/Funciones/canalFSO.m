
function [Beta0,Beta1,Beta2,Beta3,coefD,longDispersion2,longDispersion3,tauZ] = canalFSO(Dz,sigma0,h,w,c,wvl)

[Beta0,Beta1,Beta2,Beta3] = taylorBeta(h,w,c);

coefD = -2*pi*c*Beta2/(wvl*wvl);                
       

Beta0 = Beta0*1e3;                                  % [1/km] 
Beta1 = Beta1*1e15;                                 % VD [pseg/km]  
Beta2 = Beta2*1e27;                                 % GVD [pseg^2/km] 
Beta3 = Beta3*1e39;                                 % TOD [pseg^3/km] 
coefD = coefD*1e6;                                  % Coeficiente dispersión [pseg/km*nm]
longDispersion2 = 2*(sigma0^2)/abs(Beta2);          % Longitud dispersión lineal [km]  
longDispersion3 = (sqrt(2)*sigma0)^3/abs(Beta3);    % Longitud dispersión lineal [km] 
tauZ = Beta1*Dz;                                    % Retardo de grupo [pseg]

% Nos cargamos Beta0 y Beta1 porque no tienen impacto en la anchura y la
% deformación de los pulsos
Beta0 = 0;
Beta1 = 0;

