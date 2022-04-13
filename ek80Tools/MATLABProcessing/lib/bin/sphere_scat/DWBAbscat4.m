% Analytical approximation solution from echo average paper (Stanton et al, 1994)
function [ka0, ang, sigma_bs]=DWBAbscat4(para)

% shape parameters

L_a=para.shape.L_a;
L_d=L_a/2;a=para.simu.ka/para.simu.k;		% in mm
L=a*L_a;
d=2*a;

f = para.simu.f0*1e3;  % (Hz)
c = para.phy.cw;    % (m/s)

z=para.phy.g1*para.phy.hL;
r = (z-1)/(z+1);   % reflection coef. from fit of curves to data
s = 0.0;    % relative standard deviation of length

%%% compute sigmabs for each polychaete size class using tim's Bent
%%% Cylinder Model
step1 = cos(pi*f.*d.*c^(-1).*(4*(ones(size(L)))-0.5*pi*(pi*f.*d.*c^(-1)+0.4*(ones(size(L)))).^(-1)));
step2 = (1-exp(-8*pi^2*f^2.*d.^2.*s^2*c^(-2)).*step1);
sigma_bs= (0.08*r^2.*L.^2*L_d^(-1)).*step2;      % backscattering cross-section for each size class
ka0=mean(para.simu.ka);
ang=0;