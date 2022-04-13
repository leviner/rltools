%% demo program to call elastic_fs.m

clear

opt = 5;                % 1 = rigid and fixed
                        % 2 = pressure release
                        % 3 = fluid
                        % 4 = elastic
                        % 5 = elastic shell
                        
out_flag=1;				% modular of form function
proc_flag=1;			% form function vs ka
scale=1;				% linear spacing in ka
n=2000;					% number of computation points 
x0=0.1;					% starting ka value
xe=20;					% end ka value
theta=180;				% backscattering

a=30*0.025/2;
r=1;                    % ratio b/a for non-shelled objects

S=33;
T=15;
P=20;
cw0=sw_svel(S,T,P);
rho0=sw_dens(S,T,P)/1000;

T1=15;
cw1=sw_svel(S,T1,P);
rho1=sw_dens(S,T1,P)/1000;

target_index=3;
% 1: tungsten carbide
% 2: copper
% 3: aluminum
% 4: Stainless Steel
% 5: specified

switch opt
    case 1          % rigid/fixed
        para_rgd=[n x0 xe theta];
        [ka, fm]=rgd_sft_fs(proc_flag,1,scale,out_flag,para_rgd);
    case 2          % pressure release
        para_rgd=[n x0 xe theta];
        [ka, fm]=rgd_sft_fs(proc_flag,0,scale,out_flag,para_rgd);
    case 3           % fluid sphere 
        g=0.0012;              % rho2/rho1
        h=0.22;                % c2/c1
        para_fld=[n x0 xe g h theta];
        [ka, fm]=fluid_fs(proc_flag,scale,out_flag,para_fld);
    case 4           % elastic sphere
        g=14.65;				% rho2/rho1
        hc=3.7;					% cL2/c1
        hs=2.4;					% cT2/c1
        para_ela=[n x0 xe g hc hs theta];
        [ka, fm]=elastic_fs(proc_flag,scale,out_flag,para_ela);
    case 5          % elastic shelled sphere
% region 1(medium 1, r > a ); medium 2(elastic shell, b < r <= a);medium 3(inside fluid region, r <= b)

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
case 3   % Alumimum
   rho=2.7;
   cc=6260;
   cs=3080;
case 4   % Stainless Steel
   rho=7.9;cc=3.7368*cw;cs=2.08*cw;
otherwise
   rho=1.5;cc=1.5*1485;cs=0.8*1485;   
end

if 1 % water filled
        rho12=rho0/rho;         % rho1/rho2
        rho32=rho1/rho;           % rho3/rho2
        hc=cc/cw0;              % cL2/cw0
        hs=cs/cw0;              % cT2/cw0
        h3=cw1/cw0;             % c3/cw0
else  % gas filled
        rho12=rho0/rho;         % rho1/rho2
        rho32=0.00127*rho12;           % rho3/rho2
        hc=cc/cw;              % cL2/c1
        hs=cs/cw;              % cT2/c1
        h3=0.22;                % c3/c1
end
        thk=0.2*0.0254;
        r=1-thk/a;               % b/a
        para_shl=[n x0 xe rho12 rho32 hc hs h3 r theta];
        [ka, fm]=shell_fs(proc_flag,scale,out_flag,para_shl);
end
TS=20*log10(fm*a/2);
freq=1500*ka/(2*pi*a);                % frequency in kHz
figure(2)
if scale == 1
%   plot(ka,TS)
%   plot(ka,0.5*fm/sqrt(pi))
   plot(freq*1e-3,TS,'c')
else
   freq=1.5*ka*r/(2*pi);                % frequency in kHz
   semilogx(freq*1e-3,TS,'-r')
end
xlabel('ka')
xlabel('Frequency (kHz)')
ylabel('TS (dB)')
%ylabel('|f|/\sqrt{\pi}}a')
grid on
axis([0 12 -70 10])