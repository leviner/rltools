%% examples 
addpath c:\matlab\siphon\lib
addpath c:\matlab\fluidsphere

clear

addpath lib
opt=2;			% 1-DWBA,  2-Bubble,  3-Shell
para=set_para;				% obtain default parameters
tic
switch opt
case 1				% DWBA
% average over angle and length
%  euphausiids
L=[25 30 35 40];			% lengths of different animal groups in mm
para.simu.aveA_flag=1;	% average over angle
para.simu.nA=30;			% number of discrete angles to be averaged over
para.shape.ang=0;			% mean incident angle (tilt angle), 0-broadside
para.simu.aveL_flag=1;	% average over length
para.shape.Lstd=0.1;	% mean incident angle (tilt angle), 0-broadside
para.simu.nL=20;			% number of discrete angles to be averaged over
para.simu.f0=120;
para.simu.fe=120;
para.simu.df=2;
for i=1:length(L)
  para.shape.L=L(i);
  [sigma_bs_L, para]=zoo_bscat(para);
  sigma_bs(i,1:length(para.simu.freq0))=sigma_bs_L*(para.shape.L*1e-3).^2;
end
case 2
   %z=   1.5656667e+001   ;
   z=5:10:205;
   for Index=1:length(z);
para.simu.model=2;
para.simu.aveL_flag=1;	% average over length
para.simu.nL=6;			% number of discrete length to be averaged over
para.simu.f0=10;        % starting freq
para.simu.fe=500;       % ending freq
para.simu.df=1;
%L=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 1.0 1.5 2];			% bubble diameter  lengths of different animal groups in mm
   L=[ 0.21];
%%%% modified on 6-28-2001 to create para.simu.ka
a=L/2;
para.simu.ka=para.simu.k*a;
para.shape.L_a=2;
para.phy.g1=0.0012;
para.phy.g1=para.phy.g1*(1+0.1*z(Index));
para.phy.hL=0.22;
for i=1:length(L)
  para.shape.L=L(i);
  [sigma_bs_L, para]=zoo_bscat(para);
  sigma_bs(i,1:length(para.simu.freq0))=sigma_bs_L*(para.shape.L*1e-3).^2;
end
TS(Index,:)=10*log10(sigma_bs);
end 
case 3
para.simu.model=3;
para.simu.aveL_flag=1;	% average over length
para.simu.nL=100;			% number of discrete angles to be averaged over
para.simu.f0=500;
para.simu.fe=500;
para.simu.df=1;
L=[0.7 1.0 1.5 2];			% lengths of different animal groups in mm
para.shape.Lstd=0.15;
para.shape.L_a=2;
para.shape.shl=0.025;
para.phy.g1=2.646;
para.phy.hL=4.345;
para.phy.hT=1.5;
para.phy.g2=1.03;
para.phy.h2=1.02;
for i=1:length(L)
  para.shape.L=L(i);
  [sigma_bs_L, para]=zoo_bscat(para);
  sigma_bs(i,1:length(para.simu.freq0))=sigma_bs_L*(para.shape.L*1e-3).^2;
end
end
%TS=10*log10(sigma_bs);
%c1=plot(para.simu.freq0,TS,'k-')
%('FREQUENCY (kHz)')
%ylabel('TARGET STRENGTH (dB)')
%axis([min(para.simu.freq0)-1e-4 max(para.simu.freq0)+1e-4 -100 -40])
toc