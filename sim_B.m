function [dB,dF,in,out]=sim_B(R,r,th)
% Sphere without Lorentz correction. Background is free space. 


B0 = 3; % [Tesla]
% chi=0.273*4*pi; %deox blood
chi = -8e-6; % susceptibility of water
mu0 = 1.2566e-6; % permeability of free space

in=0; out=0;

if r<=R % inside
    dB=(2*chi/(3+chi))*B0;
    M=((3*chi)/(3+chi))*(B0/mu0);   % Lorentz correction
    dB=dB-(((2*mu0)/3)*M);
    %in=in+1;
else % outside
    dB=(chi/(3+chi))*((R/r)^3)*(3*(cos(th))^2-1)*B0;
    % No Lorentz correction - susc. of free space = 0
    %out=out+1;
end

gamma = 42.57e6; % [Hz/T]
dF=gamma*dB;