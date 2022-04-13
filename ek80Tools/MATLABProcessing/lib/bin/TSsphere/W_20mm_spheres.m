function [f,TS] = W_20mm_spheres(theta);

%% demo program to call fluid_fs.m
%% 

addpath functions
%addpath sw_property

%% type: help elastic_fs in the Matlab command window to get help information
%%       for function  elastic_fs.m

target_index=1;
% 1: tungsten carbide
% 2: copper
% 3: aluminum
proc_flag=1;			% output vs ka
scale=1;					% linear spacing in ka
out_flag=1;				% modular of form function
n=2000;					% number of computation points 
ka0=0.1;		        % smallest ka
ka1=40;					% largest ka

a=0.5*[20]*1e-3;       % radius in m (put diameter in mm in the square parenthesis)
freq=[38 120 200]*1e3;  % three frequencies in kHz ((put freq. in kHz in the square parenthesis)

cw_opt=2;               % how to get sound speed in seawater 
                        % cw_opt = 1: provide sound speed 
                        %        = 2: compute sound speed based on Salinity
                        %        (ppt), Temperature (deg C), and Pressure (dbar)
                       
if cw_opt == 1
    cw=1500;            %
    rhow=1.025;
else
    S=12;
    T=10;
    P=1;
    cw=sw_svel(S,T,P);     
    rhow=sw_dens(S,T,P)/1000;
end
switch target_index
case 1
%% properties of Tungsten carbide - CRC handbook of Chemistry and Physics, David R. Lite, editor in chief
%%                                - 77th edition, 1996-1997, 14-36
 % tungsten carbide   rho=14.9  cc=6853  cs=4171
   rho=14.9;				% density 
   cc=6853;					% speed of compressional wave 
   cs=4171;					% speed of shear wave 
case 2  % copper  rho=8.947  cc=4760   cs=2288.5
   rho=8.947;
   cc=4760;
   cs=2288.5;
case 3
   rho=2.7;
   cc=6260;
   cs=3080;
end
%theta=180;				% scattering angle in degrees (180 for backscattering)


g=rho/rhow;
hc=cc/cw;
hs=cs/cw;

ka0_t1=2*pi*freq*a(1)/cw;
indx1=find( ka0_t1 > ka0 & ka0_t1 < ka1);
ka_t1=ka0_t1(indx1);

   
para_elastic=[n ka0 ka1 g hc hs theta];
[ka, fm]=elastic_fs(proc_flag,scale,out_flag,para_elastic);
RTS=20*log10(abs(fm));
RTS_t1=interp1(ka,RTS,ka_t1);

TS=RTS+20*log10(a/2);
f=ka*cw./2./pi./a;

plot(f/1000,TS)