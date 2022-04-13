function [TS,f]=WC_TS(D,plot_flag);

%pwd_home=pwd;
%addpath([pwd_home '\bin'])
%addpath([pwd_home '\bin\seawater'])
%addpath([pwd_home '\bin\sphere_scat'])
%clear pwd_home

%fprintf('Assumes the diameter is in mm\n')
ka_min=0.05;
ka_max=40;

R=D*1e-3./2;% Sphere Radius
    % % % fprintf('Radius %2.2f mm\n',R*1000)
proc_flag=1;%ka
scale_flag=1;%linear
out_flag=3;%scaled scat amplitude
cw=sw_svel(0,20.4,3);%cw=1490;
g=14900./1000;% Density of Sphere

hL=6853/cw;
hL=6848./cw; %Longitudinal Sound Speed of Sphere
%hL=6875./cw; %Longitudinal Sound Speed of Sphere

hT=4171/cw;
hT=4161./cw; %Transverse Sound Speed of Sphere
%hT=4177./cw;

%hL=6875./cw; %Longitudinal Sound Speed of Sphere
%hT=4150./cw; %Transverse Sound Speed of Sphere
 
 % Gavin McCauley: hL = 6853 m/s, hT = 4171 m/s, and density = 14900 kg/m^3

para_flag=[10000,ka_min,ka_max,g,hL,hT,180];
[outx,outy]=elastic_fs(proc_flag,scale_flag,out_flag,para_flag);
f=cw.*outx./2./pi./R;
TS=10.*log10(outy./2.*R.^2);

if plot_flag==1
%   figure(1)
%   plot(outx,10.*log10(outy./pi./2),'k')
%   set(gca,'linewidth',[2],'fontsize',[12])
%   xlabel('ka','fontsize',[12])
%   ylabel('REDUCED TARGET STRENGTH (adB)','fontsize',[12])
%   ylim([-40 0])
figure(2)
plot(f/1000,TS,'k','linewidth',2)
set(gca,'linewidth',[2],'fontsize',[12])
xlabel('FREQUENCY (kHz)','fontsize',[12])
ylabel('TARGET STRENGTH (dB)','fontsize',[12])
title('Tungsten Carbide','fontsize',[12])
ylim([-70 -30])
grid on
end

