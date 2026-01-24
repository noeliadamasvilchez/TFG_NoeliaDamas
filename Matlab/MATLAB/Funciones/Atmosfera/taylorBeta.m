

function [Beta0,Beta1,Beta2,Beta3] = taylorBeta(h,w0,c)

% h altura del enlace en metros
% w0 frecuencia

% Constante de propagación en función de omega
betaFunction = @(w) (w/c)*refractiveIndex(2*pi*c/w,h);
%betaFunction = @(w) 3 + 4*w + 5*w*w;

% Parámetros Beta0
Beta0 = betaFunction(w0);

% Parámetro Beta1
syms w
betaFunction1 = diff(betaFunction(w),w,1);
subs(betaFunction1,w,w0);
Beta1 = double(ans);

% Parámetro Beta2
syms w
betaFunction2 = diff(betaFunction(w),w,2);
subs(betaFunction2,w,w0);
Beta2 = double(ans);

% Parámetro Beta3
syms w
betaFunction3 = diff(betaFunction(w),w,3);
subs(betaFunction3,w,w0);
Beta3 = double(ans);