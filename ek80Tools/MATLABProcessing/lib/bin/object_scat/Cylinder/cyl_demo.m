clear

opt = 3;                % 1 = rigid and fixed
                        % 2 = pressure release
                        % 3 = fluid
                        % 4 = elastic
                        % 5 = elastic shell
                        
out_flag=1;				% modular of form function
proc_flag=1;			% form function vs ka
scale=2;				% linear spacing in ka
n=1000;					% number of computation points 
x0=-3;					% starting ka value
xe=1.2;					% end ka value
theta=180;				% backscattering

a=10e-3;
r=1;                    % ratio b/a for non-shelled objects
switch opt
    case 1          % rigid/fixed
        para_rgd=[n x0 xe theta];
        [ka, fm]=rgd_sft_fc(proc_flag,1,scale,out_flag,para_rgd);
    case 2          % pressure release
        para_rgd=[n x0 xe theta];
        [ka, fm]=rgd_sft_fc(proc_flag,2,scale,out_flag,para_rgd);
    case 3           % fluid cylinder 
        g=0.0012;              % rho2/rho1
        h=0.22;                % c2/c1
        para_fld=[n x0 xe g h theta];
        [ka, fm]=fluid_fc(proc_flag,scale,out_flag,para_fld);
    case 4           % elastic cylinder
        g=14.65;				% rho2/rho1
        hc=3.7;					% cL2/c1
        hs=2.4;					% cT2/c1
        para_ela=[n x0 xe g hc hs theta];
        [ka, fm]=elastic_fs(proc_flag,scale,out_flag,para_ela);
    case 5          % elastic shelled cylinder
% region 1(medium 1, r > a ); medium 2(elastic shell, b < r <= a);medium 3(inside fluid region, r <= b)
        rho12=1/1.026;          % rho2/rho1
        rho32=0.0012/rho12;     % rho3/rho1
        hc=1.0;               % cL2/c1
        hs=1e-16;           % cT2/c1
        h3=0.22;                % c3/c1
        r=0.99999;               % b/a
        para_shl=[n x0 xe rho12 rho32 hc hs h3 r theta];
        [ka, fm]=shell_fc(proc_flag,scale,out_flag,para_shl);
end
TS=20*log10(fm*2*a);
figure(1)
if scale == 1
   plot(ka,TS)
else
   freq=1.5*ka*r/(2*pi);                % frequency in kHz
   semilogx(freq,TS,'-r')
end
xlabel('ka')
ylabel('TS (dB)')
grid
