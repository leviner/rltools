% % Stanton's high freq. Kirchhoff approximation
function [ka0, ang, sigma_bs]=DWBAbscat5(para)


L_a=para.shape.L_a;
L_d=L_a/2;                  % ratio of length to width
a=para.simu.ka/para.simu.k;		% in mm
L=a*L_a;                    % animal length
d=L;                        % diameter of medusa bell, or width of the animal
z=para.phy.g1*para.phy.hL;  % acoustic impedance
R=(z-1)/(z+1);              % plane wave refelction coefficient

% compute sigmabs for each individual gelatinous plankter using Gelatinous model
% step one--calculate principal radii of curvature
r1=d.*.25;
r2=d.*.75;	

%step two -- compute high frequncy backscattering cross section based on
%the Kirchhoff approximation
sigma_bs=.25*((r1.^2)+(r2.^2))*(R.^2);  % backscattering cross-section for each individual
ka0=mean(para.simu.ka);
ang=0;