function [TS,f]=SS_TS(D,plot_flag);

addpath c:\matlab\TSsphere\functions
addpath c:\matlab\object_scat\Sphere
addpath c:\matlab\toolboxes\seawater
addpath c:\matlab\sphere_scat\

%fprintf('Assumes the diameter is in mm\n')
ka_min=0.01;
ka_max=50;

R=D*1e-3./2;% Sphere Radius
fprintf('Radius %2.2f mm\n',R*1000)
proc_flag=1;%ka
scale_flag=1;%linear
out_flag=3;%scaled scat amplitude

cw=sw_svel(0,20.4,3);%cw=1490;
g=7959./1000;% Density of Stainless Steel Sphere +/- 1 km/m3
hL=5669.94/cw; %Longitudinal/Compressional Sound Speed of Sphere +/-21.71m/s
hT=2987.69/cw; %Transverse/Shear Sound Speed of Sphere +/- 2.24 m/s


para_flag=[10000,ka_min,ka_max,g,hL,hT,180];
[outx,outy]=elastic_fs(proc_flag,scale_flag,out_flag,para_flag);
f=cw.*outx./2./pi./R;
TS=10.*log10(outy./2.*R.^2);

if plot_flag==1;
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