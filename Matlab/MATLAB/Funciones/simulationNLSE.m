

function [E_out,F_out,fourierF_out] = simulationNLSE(F_in,Dz,alpha,Beta1,Beta2,Beta3,Carrier_out,omega,n,nn)

% The evolution of slow varying complex envelopes F(z,t) of optical pulses along a single-mode optical
% fiber is governed by the nonlinear Schrödinger equation (NLSE)

F_out_sim = zeros(nn,n);        % Propagation matrix 

% Pulso en z = 0
F_out_sim(1,:) = F_in;          % F0*F_IN(t)
h = Dz/nn;                      % Spatial step (propagation step) 

% Dispersion operator in Frequency space: Ecuación Shordinger No Lineal 
D = -(-1i*0.5*Beta2*(1i*omega).^2 + (1/6)*Beta3*(1i*omega).^3 - alpha/2); % Signo menos lo pongo yo

% Split-Step Fourier Algorithim 
for ii = 2:nn  
    %N = -1i*gamm*abs(F_out_sim(ii-1,:)).^2;  % Signo menos lo pongo yo
    %F_out_sim(ii,:) = exp(h*N).*ifft(fftshift(exp(h*D).*fftshift(fft(F_out(ii-1,:)))));
    F_out_sim(ii,:) = ifft(fftshift(exp(h*D).*fftshift(fft(F_out_sim(ii-1,:)))));
end

% Desplazamos para incorporar retardo de grupo: Beta1*Dz
if Beta1 == 0
    F_out = F_out_sim(end,:);               % Pulso en z = Dz
else
    ms = floor(Beta1*Dz/dt);
    F_out = circshift(F_out_sim(end,:),ms); % Pulso en z = Dz
end
E_out = F_out.*Carrier_out;

% Fourier Space of the propagation steps
fourierF_out = zeros(size(F_out_sim)); 
for ii = 1:nn
    fourierF_out(ii,:) = fftshift(fft(F_out_sim(ii,:)));
end

