function [ka0, f]=Bubblescat1(para)
%% compute backscattering by a gas bubble

proc_flag=1;
scale=1;
out_flag=1;
x0=min(para.simu.ka1);
xe=max(para.simu.ka1);
ns=length(para.simu.ka1);
g=para.phy.g1;h=para.phy.hL;
model_para=[ns x0 xe g h 180];		% backscattering
[ka0, f0]=fluid_fs(proc_flag,scale,out_flag,model_para);
f=f0/4;		% f*L=f*2a=a/2*f0 -> scattering amplitude

