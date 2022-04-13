%% demo program to call elastic_fs.m

clear

out_flag=1;				% modular of form function
proc_flag=1;			% form function vs ka
scale=1;					% linear spacing in ka
n=1000;					% number of computation points 
x0=0.1;					% starting ka value
xe=10;					% end ka value
g=2.65;					% density ratio
hc=3.7;					% compressional sound speed contrast
hs=2.4;					% sheer speed contrast
theta=180;				% backscattering

para=[n x0 xe g hc hs theta];

[ka, fm]=elastic_fs(proc_flag,scale,out_flag,para);

plot(ka,fm)
xlabel('ka')
ylabel('Form Function |f|')
