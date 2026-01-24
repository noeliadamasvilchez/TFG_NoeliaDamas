
function out = TPS_OTOPS(args)

% % Script que calcula los espectros individuales OTOPS para temperatura,
% salinidad, co-espectro y espectro total. Uso de TEOS-10 toolbox.

% ENTRADAS (contenidas en args)
%    SP o Salinidad Práctica [PSU], T [°C], Presión [dbar], Longitud de onda [nm],
%    ε o tasa de disipación de energía cinética [10^-4 o 1e-4 m^2/s^3], beta0 o cte Obukhov-Corrsin (0.72), 
%    Xt [10^-6 o 1e-6 K^2/s] y H_grad [-20 °C/ppt].
%    
% Opcionales
%   lat, lon                 % por defecto 0,0 para evitar NaN en TEOS
%   alpha, beta              % extraer de TEOS


% SALIDA (contenida en el struct out)
%   Coeficientes A y B, SA o Salinidad Absoluta, μ o viscosidad dinámica,
%   ρ o densidad, η o escala de Kolmogorov, σT o conductividad térmica,
%   cP o calor específico, Pr o número de Prandtl, Sc o número de Schmidt, 
%   cT, cS y cTS o coeficientes adimensionales, Xt, Xs, XTS o tasas de
%   disipación y Phi_T, Phi_S, Phi_TS o funciones que devuelven los espectros parciales.


% Chequeo de las entradas (asegurar que están todos los argumentos)
req = {'SP','T','lambda_nm','epsilon','beta0','Xt','H_grad'};
for i=1:numel(req)
    if ~isfield(args,req{i}), error('Falta args.%s',req{i}); end
end

% Lectura de entradas
SP          = args.SP;
T           = args.T;
p           = args.Pdbar;
lambda_nm   = args.lambda_nm;
epsilon     = args.epsilon;
beta0       = args.beta0;
Xt          = args.Xt;     
H_grad      = args.H_grad;
lat         = args.lat;
lon         = args.lon;

% 1) Valores a obtener de TEOS-10: SA, CT o Temperatura Conservativa, 
% alpha, beta 

SA = gsw_SA_from_SP(SP, p, lon, lat);       % [g/kg] ≈ [ppt]
CT = gsw_CT_from_t(SA, T, p);               % [°C]
alpha = gsw_alpha(SA, CT, p);               % [1/K]
beta  = gsw_beta(SA, CT, p);                % [kg/g]
   
out.SA    = SA;
out.CT    = CT;
out.alpha = alpha;
out.beta  = beta;


% 2) Coeficientes A y B 

A = (-1.05e-6*SA) + (3.2e-8*CT*SA) - (4.04e-6*CT) - (0.00423./lambda_nm);   % [K^-1]
B =  1.779e-4 - 1.05e-6*T + 1.6e-8*T.^2 + 0.01155./lambda_nm;               % [ppt^-1]

out.A = A;
out.B = B;

% 3) Viscosidad μ y Densidad ρ para escala de Kolmogorov

% Cálculo del primer término, F1
F1 = [0.15700386464 * (CT + 64.992620050)^2 - 91.296496657]^-1 + 4.2844324477e-5;
F1_denominador = 0.15700386464 * (CT + 64.992620050)^2 - 91.296496657;
F1 = (F1_denominador)^(-1) + 4.2844324477e-5;

% Cálculo del segundo término, F2
F2_factor_S = 1.5409136040e-3 + 1.9981117208e-5*CT - 9.5203865864e-8*(CT^2);
F2_factor_S2 = 7.9739318223e-6 - 7.5614568881e-8*CT + 4.7237011074e-10*(CT^2);
F2 = 1 + F2_factor_S * SA + F2_factor_S2 * (SA^2);

% Obtención de la viscosidad
mu = F1 * F2;
out.mu = mu;


% Parte inicial de la fórmula que solo depende de T.
PT = 9.9992293295e+2 ...
     + 2.0341179217e-2 * CT ...
     - 6.1624591598e-3 * (CT^2) ...
     + 2.2614664708e-5 * (CT^3) ...
     - 4.6570659168e-8 * (CT^4);

% Parte de la fórmula que multiplica a S (y también depende de T).
PS_factor_interno = 8.0200240891e-1 ...
                    - 2.0005183488e-3 * CT ...
                    + 1.6771024982e-5 * (CT^2) ...
                    - 3.0600536746e-8 * (CT^3);

PS = SA * PS_factor_interno ...
     - 1.6132224742e-11 * (CT^2) * SA; 

% Obtención de la densidad
rho = PT + PS;
out.rho = rho;

% 4) Escala de Kolmogorov η
eta = mu.^(3/4) .* rho.^(-3/4) .* epsilon.^(-1/4);   % [m]
%kef = 1./eta; 
out.eta = eta;
%out.kef = kef;

% 5) Conductividad térmica σT y Calor específico cP 

% Temperatura en Kelvin
TK = T + 273.15;

% Coeficientes de salinidad
T0 = 5.328 - 9.76e-2 * SA + 4.04e-4 * (SA^2);
T1 = (-6.913e-3 + 7.351e-4 * SA - 3.15e-6 * (SA^2));
T2 = (9.6e-6 - 1.927e-6 * SA + 8.23e-9 * (SA^2));
T3 = (2.5e-9 + 1.666e-9 * SA - 7.125e-12 * (SA^2));

% Expresión corchete
Corchete = T0 + T1 * TK + T2 * (TK^2) + T3 * (TK^3);

% Obtención del calor específico
cp = 1000 * Corchete; 
out.cp = cp;

% Obtención conductividad térmica
Th = 1.00024 .* CT;
Sh = SA ./ 1.00472;
Tk = Th + 273.15;
log10_sigma = log10(240 + 0.0002.*Sh) ...
            - 3 ...
            + 0.434 .* ( 2.3 - (343.5 + 0.037.*Sh)./Tk ) ...
                    .* ( 1 - Tk./(647.3 + 0.03.*Sh) ).^(1/3);

sigma_T = 10.^log10_sigma;          % [W/(m·K)]
out.sigma_T = sigma_T;

% 6) Pr, Sc y c_i 

Pr  = cp.*mu ./ sigma_T;
Sc  = (mu.^2) ./ (5.954e-15 * CT .* rho);
out.Pr = Pr;
out.Sc = Sc;

cT  = (0.072)^(4/3) * beta0 ./ Pr;
cS  = (0.072)^(4/3) * beta0 ./ Sc;
cTS = (0.072)^(4/3) * beta0 .* (Pr + Sc) ./ (2*Pr.*Sc);
out.cT = cT;
out.cS = cS;
out.cTS = cTS;

% 7) Xt, Xs, XTS 

 Rp = abs(H_grad)*alpha*(beta^-1);
 out.Rp = Rp;

 dr = Rp + (Rp^0.5) * ((Rp - 1)^0.5);
 out.dr = dr;

 Xs  = dr * Xt /  H_grad^2 ;            % [ppt^2/s]
 out.Xs = Xs;

 XTS = 0.5 * (1 + dr) * Xt / H_grad;    % [K·ppt/s] 
 out.XTS = XTS;


% 8) Espectros parciales 

prefijo = beta0 * epsilon.^(-1/3) / (4*pi);

Phi_T_fun  = @(k) prefijo .* ...
    ( 1 + 21.61.*(k.*eta).^0.61 .* cT.^0.02 ...
        - 18.18.*(k.*eta).^0.55 .* cT.^0.04 ) ...
    .* k.^(-11/3) .* Xt ...
    .* exp( -174.90 .* (k.*eta).^2 .* cT.^0.96 );

Phi_S_fun  = @(k) prefijo .* ...
    ( 1 + 21.61.*(k.*eta).^0.61 .* cS.^0.02 ...
        - 18.18.*(k.*eta).^0.55 .* cS.^0.04 ) ...
    .* k.^(-11/3) .* Xs ...
    .* exp( -174.90 .* (k.*eta).^2 .* cS.^0.96 );

Phi_TS_fun = @(k) prefijo .* ...
    ( 1 + 21.61.*(k.*eta).^0.61 .* cTS.^0.02 ...
        - 18.18.*(k.*eta).^0.55 .* cTS.^0.04 ) ...
    .* k.^(-11/3) .* XTS ...
    .* exp( -174.90 .* (k.*eta).^2 .* cTS.^0.96 );

Phi_n_fun  = @(k) (A.^2).*Phi_T_fun(k) + (B.^2).*Phi_S_fun(k) + 2*A.*B.*Phi_TS_fun(k);

% Salida de funciones:
out.Phi_T_fun  = Phi_T_fun;     % Espectro de temperatura
out.Phi_S_fun  = Phi_S_fun;     % Espectro de salinidad
out.Phi_TS_fun = Phi_TS_fun;    % Espectro del co-espectro
out.Phi_n_fun  = Phi_n_fun;     % Espectro total

end
