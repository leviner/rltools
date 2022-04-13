function [ka0, f]=Shellscat1(para)

proc_flag=1;
scale=1;
out_flag=1;
x0=min(para.simu.ka1);
xe=max(para.simu.ka1);
ns=length(para.simu.ka1);
g12=1/para.phy.g1;
g31=para.phy.g2;
g32=g31*g12;
hL=para.phy.hL;
hT=para.phy.hT;
h2=para.phy.h2;
r=1-para.shape.shl;
model_para=[ns x0 xe g12 g32 hL hT h2 r 180];		% backscattering
[ka0, f0]=shell_fs(proc_flag,scale,out_flag,model_para);
f=f0/4;		% f*L=f*2a=a/2*f0 -> scattering amplitude
