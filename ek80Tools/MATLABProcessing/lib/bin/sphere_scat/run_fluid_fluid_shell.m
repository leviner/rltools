clear
close all
clc

f1=1e3;
f2=1000e3;

c1=1476;
c3=1478;
c2=1540;

rho1=993;
rho3=997;
rho2=1045;

h21=c2/c1;
h31=c3/c1;
g21=rho2/rho1;
g31=rho3/rho1;

a=0.02e-3;
b=5e-3;

kb1=2*pi*f1./c1.*b;
kb2=2*pi*f2./c1.*b;

N=500;

para1=[N,kb1,kb2,rho1,rho2,rho3,c1,c2,c3,a,b,180];

%% CYLINDRICAL SHELL
[kb_C,fbs_shell_C]=fluid_fluid_shell_cyl_new(1,1,3,para1);
sigma_shell_C=(fbs_shell_C./pi).^2;

semilogx(kb_C,10*log10(sigma_shell_C),'k','linewidth',[2])

g=rho2./rho1;
h=c2./c1;
temp=(1-g.*h.^2)./(2.*g.*h.^2) + (1-g)./(1+g);
temp1=0.25*kb_C.^4.*temp.^2;
hold on
semilogx(kb_C,10*log10(temp1))
hold off

%% CYLINDER
[kb_C,fbs_shell_C]=fluid_cyl(1,1,3,para1);
sigma_shell_C=(fbs_shell_C./pi).^2;

hold on
semilogx(kb_C,10*log10(sigma_shell_C),'r','linewidth',[1])
hold off
axis([0.01 80 -70 10])

return

%%SPHERICAL SHELL
[kb,fbs_shell]=fluid_fluid_shell(1,1,3,para1);
fbs_shell=fbs_shell.*b;
return
%% SPHERE B MATERIAL PROPERTIES 2
para2=[N,kb1,kb2,rho2./rho1,c2./c1,180];
[kb_sphr_b,fbs_sphr]=fluid_fs(1,1,3,para2);
fbs_sphr_b=fbs_sphr.*sqrt(pi).*b;

%% SPHERE RADIUS A MATERIAL PROPERTIES 3
para3=[N,kb1,kb2,rho3./rho1,c3./c1,180];
[kb_sphr_a,fbs_sphr]=fluid_fs(1,1,3,para3);
fbs_sphr_a=fbs_sphr.*sqrt(pi).*a;


plot(kb,10*log10(fbs_shell.^2),'k','linewidth',[2])
hold on
%plot(kb_C,10*log10(fbs_shell_C.^2),'g','linewidth',[2])
plot(kb_sphr_b,10*log10(fbs_sphr_b.^2),'r','linewidth',[2])
plot(kb_sphr_a,10*log10(fbs_sphr_a.^2),'b','linewidth',[2])
hold off
ylim([-130 -50])
xlabel('k*radius','fontsize',[12])
ylabel('TS (dB)','fontsize',[12])
set(gca,'linewidth',[2],'fontsize',[12])


return
%legend('Fluid-filled Fluid Shell (Radii: inner=a, outer=b)',...
 %   'Fluid Sphere Radius b(>a) Material Properties II',...
 %   'Sphere Radius a Material Properties III',1)

 
 
 
 
 
 
 
%%% REGULAR DWBA SPHERE RADIUS B Properties 2
kappa2=1./(rho2.*c2.^2);
kappa1=1./(rho1.*c1.^2);

gamma_kappa21=(kappa2-kappa1)./kappa1;
gamma_rho21=(rho2-rho1)./rho2;

k1=kb_sphr_b./b;
k2=k1.*c1./c2;
temp1=-2*k1.*b.*cos(2*k2*b)+sin(2*k2*b);
fbs_DWBA_b_II=(gamma_kappa21-gamma_rho21)./8./(k1).*temp1;

%%% REGULAR DWBA SPHERE RADIUS A Properties 3
kappa3=1./(rho3.*c3.^2);
kappa1=1./(rho1.*c1.^2);

gamma_kappa31=(kappa3-kappa1)./kappa1;
gamma_rho31=(rho3-rho1)./rho3;

k3=k1.*c1./c3;
temp2=-2*k1.*a.*cos(2*k3*a)+sin(2*k3*a);
fbs_DWBA_a_III=(gamma_kappa31-gamma_rho31)./8./k1.*temp2;

%%% DWBA SHELL ONLY of Material Properties II
fbs_DWBA_shell_II=(gamma_kappa21-gamma_rho21)./8./(k1).*(temp1-temp2);


hold on
plot(kb_sphr_b,10*log10(fbs_DWBA_b_II.^2),'m')
plot(kb_sphr_b.*a/b,10*log10(fbs_DWBA_a_III.^2),'c')
plot(kb_sphr_b,10*log10(fbs_DWBA_shell_II.^2),'g')
hold off