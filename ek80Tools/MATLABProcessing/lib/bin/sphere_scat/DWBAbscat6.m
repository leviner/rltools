function [ka0, ang, sigma_bs]=DWBAbscat6(para)
% Spherical Elastic Shell Model (Pteropod) from Stanton et al 1994

L_a=para.shape.L_a;
L_b=para.shape.L_b;
L_d=L_a/2;
a=para.simu.ka/para.simu.k;		% in mm
d=a*sqrt(L_a*L_b);              % (m) average radius of individual pteropods in 10 samples (nans fill out short columns)

f = para.simu.f0*1000;  % (Hz)
c = para.phy.cw;    % (m/s)
z=para.phy.g1*para.phy.hL;  % acoustic impedance
R=(z-1)/(z+1);              % plane wave refelction coefficient

% compute sigmabs for each individual pteropod using tim's Elastic Shell model

step0=d.^4;
step1 = (1+(25/10)*pi^4*f^4*c^(-4).*step0).^(-1);
step2 = ((25/144)*pi^4.*d.^6*f^4*R^2*c^(-4)).*step1;
sigma_bs=step2;         % backscattering cross-section for each individual pteropod

ka0=mean(para.simu.ka);
ang=0;