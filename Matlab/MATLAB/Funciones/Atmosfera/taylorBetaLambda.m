

function [Beta0,Beta1,Beta2,Beta3] = taylorBetaLambda(h,wvl0)

% h altura del enlace en metros
% wvl0 longitud de onda

% Constante de propagación en función de lambda
betaFunction = @(wvl) (2*pi/wvl)*refractiveIndex(wvl,h);

% Parámetros Beta0
Beta0 = betaFunction(wvl0);

% Parámetro Beta1
syms wvl
betaFunction1 = diff(betaFunction(wvl),wvl,1);
subs(betaFunction1,wvl,wvl0);
Beta1 = double(ans);

% Parámetro Beta2
syms wvl
betaFunction2 = diff(betaFunction(wvl),wvl,2);
subs(betaFunction2,wvl,wvl0);
Beta2 = double(ans);

% Parámetro Beta3
syms wvl
betaFunction3 = diff(betaFunction(wvl),wvl,3);
subs(betaFunction3,wvl,wvl0);
Beta3 = double(ans);