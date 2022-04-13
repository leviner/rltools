function		para=set_para


%% geometric parameters
para.shape.L=.21;					% mean length in mm
para.shape.Lstd=0.1;				% relative standard deviation std(L)/L
para.shape.L_a=16;				% L/a
para.shape.ang=0;					% mean incident angle in degree
para.shape.dang=20;				% standard deviation std(ang)
para.shape.rho_L=3;				% rho/L
para.shape.shl=0;					% relative shell thickness 
para.shape.dshl=0;				% relative standard deviation in shell thickness
para.shape.order=10;				% parameter controls tapering funtion

%% physical parameters
para.phy.g1=1.0357;				% density contrast: rho1/rhow (outer layer)
para.phy.g2=1.0357;				% density contrast: rho2/rhow (inner layer)
para.phy.hL=1.0279;				% compressional sound speed contrast (outer layer)
para.phy.hT=0;						% transverse sound speed contrast (outer layer)
para.phy.h2=1.0279;				% compressional sound speed contrast (inner layer)
para.phy.cw=1500;					% sound speed in water (m/s)

%% simulation parameters
para.simu.f0=120;						% start frequecy in kHz
para.simu.fe=120;						% end frequency in kHz
para.simu.df=1;						% frequency increment in khz
para.simu.model=1;					% scatt. model: 1-fluid, 2-gas, 3-shell
para.simu.aveL_flag=0;			% average over length flag
para.simu.aveL_PDF=2;				% PDF model: 1-uniform, 2-Gaussian, 
para.simu.nL=1;						% number of discrete L's to be averaged
para.simu.aveA_flag=0;			% average over angle flag
para.simu.aveA_PDF=1;				% PDF model: 1-uniform, 2-Gaussian, 
para.simu.nA=1;						% number of discrete A's to be averaged
para.simu.freq0=para.simu.f0:para.simu.df:para.simu.fe;
para.simu.freq1=para.simu.freq0;
para.simu.min_ni=50;				% minimum integration points
if para.simu.aveA_flag == 1
	para.simu.ang=linspace(para.shape.ang-3*para.shape.dang, ...
							  para.shape.ang+3*para.shape.dang,para.simu.nA);
else
	para.simu.ang=para.shape.ang;
end

ns=para.simu.min_ni;				% number of sample points per wave length
para.simu.k=2*pi*para.simu.freq1*1e3/para.phy.cw;		% wave number
kLmax=max(para.simu.k)*para.shape.L*1e-3*1.3;		% 1+ 3 std(dL/L)
para.simu.ni=max(para.simu.min_ni,ceil(kLmax*ns/(2*pi)));		% integration points along z-axis

